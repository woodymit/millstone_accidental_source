"""
Methods for filtering over Variants.

The fields that we allow querying against vary across different models,
and for each model, the field may either be a column in the SQL database, a key
in the catch-all serialized key-value dictionary.

A naive implementation would pull all Variant objects from the SQl database
(limited to the ReferenceGenome provided) into memory, and then iterate over
all these objects for each condition in the query. However, we could reduce
the number of objects pulled into memory by extracting the parts of the query
that correspond to columns in the SQL table. The strategy for doing this is
to convert the query into disjunctive normal form, an OR of ANDs, and then
evaluate each conjunction (AND clause) separately. Finally, take the union
of the result of each evaluated conjunction.

We can use sympy, a python symbolic mathematics library, to help convert the
user's query into DNF. The main idea here is to take the raw query, replace
conditions with a symbol, convert to DNF using sympy, and then evaluate each
conjunctive clause as described above.

Other implementation nuances to note:
    * Since we are comparing across multiple tables that are inter-related,
      remember to use Django's select_related() method where appropriate, in
      order to limit the number of independent hits of the SQL db.
"""

from collections import defaultdict
from django.db.models import Q
from sympy.logic import boolalg

from main.models import ExperimentSample
from main.models import Region
from main.models import VariantAlternate
from main.models import VariantEvidence
from variants.common import ALL_SQL_KEY_MAP_LIST
from variants.common import DELIM_TO_Q_POSTFIX
from variants.common import EXPRESSION_REGEX
from variants.common import SAMPLE_SCOPE_REGEX
from variants.common import SAMPLE_SCOPE_REGEX_NAMED
from variants.common import GENE_REGEX
from variants.common import GENE_REGEX_NAMED
from variants.common import SET_REGEX
from variants.common import SET_REGEX_NAMED
from variants.common import TYPE_TO_SUPPORTED_OPERATIONS
from variants.common import get_all_key_map
from variants.common import get_delim_key_value_triple
from variants.common import ParseError



###############################################################################
# Helper objects.
###############################################################################

def symbol_generator():
    """Generator that yields characters in alphabetical order, starting with
    'A' and continuing on to 'Z', followed by 'a' and ending with 'z'.

    NOTE: The current implementation runs out of symbols after 26 * 2 symbols.
    """
    final_symbol_ord = ord('z')
    current_char = 'A'
    while True:
        yield current_char
        next_ord = ord(current_char) + 1
        if next_ord > final_symbol_ord:
            raise StopIteration()
        current_char = chr(next_ord)


class FilterEvalResult(object):
    """Wraps the result of evaluating a filter condition.

    Provides utility methods for combining results.
    """

    def __init__(self, variant_set, variant_id_to_metadata_dict):
        if not isinstance(variant_set, set):
            variant_set = set(variant_set)
        self.variant_set = variant_set
        self.variant_id_to_metadata_dict = variant_id_to_metadata_dict

    def __or__(self, other):
        return self.combine(other, '|')

    def __and__(self, other):
        return self.combine(other, '&')

    def combine(self, other, op_string):
        """Method that returns a new FilterEvalResult that is the combination
        of this one and other.

        Args:
            other: The FilterEvalResult to combine with.
            op_string: Either '&' or '|'.

        Returns:
            A new FilterEvalResult object.
        """
        assert isinstance(other, FilterEvalResult)
        assert op_string in ['&', '|']

        # Combine the Variant sets.
        if op_string == '&':
            new_variant_set = self.variant_set & other.variant_set
        elif op_string == '|':
            new_variant_set = self.variant_set | other.variant_set
        else:
            raise AssertionError("Unsupported op: %s" % op_string)

        # Build up the new metadata map.
        new_variant_id_to_metadata_dict = {}
        for variant in new_variant_set:
            merged_filter_metadata = {}

            self_filter_metadata = self.variant_id_to_metadata_dict.get(
                    variant.id, {})
            other_filter_metadata = other.variant_id_to_metadata_dict.get(
                    variant.id, {})

            # Merge passing sample ids.
            self_passing_genomes = self_filter_metadata.get(
                    'passing_sample_ids', set())
            other_passing_genomes = other_filter_metadata.get(
                    'passing_sample_ids', set())
            if op_string == '&':
                merged_filter_metadata['passing_sample_ids'] = (
                        self_passing_genomes & other_passing_genomes)
            else:
                merged_filter_metadata['passing_sample_ids'] = (
                        self_passing_genomes | other_passing_genomes)

            # Save the filter metadata.
            new_variant_id_to_metadata_dict[variant.id] = merged_filter_metadata

        return FilterEvalResult(new_variant_set,
                new_variant_id_to_metadata_dict)


FILTER_SCOPE__ALL = 'ALL'
FILTER_SCOPE__ANY = 'ANY'
FILTER_SCOPE__ONLY = 'ONLY'
VALID_FILTER_SCOPES = set([
    FILTER_SCOPE__ALL,
    FILTER_SCOPE__ANY,
    FILTER_SCOPE__ONLY
])

class FilterScope(object):
    """Represents the scope that a filter should be applied over.
    """

    def __init__(self, scope_type, sample_ids):
        """
        Args:
            sample_ids: Set of sample ids.
            scope_type: A scope in VALID_FILTER_SCOPES.
        """
        assert scope_type in VALID_FILTER_SCOPES, "Invalid scope type."

        self.sample_id_set = set(sample_ids)
        self.scope_type = scope_type


    @classmethod
    def parse_sample_ids(clazz, sample_id_string):
        """Turns a comma-separated list of sample uids or names into ids.
        """
        sample_uids_or_names = sample_id_string.split(',')
        sample_uids_or_names = [s.strip() for s in sample_uids_or_names]
        sample_ids = [ExperimentSample.objects.get(uid=uid).id for uid
                in sample_uids_or_names]
        return sample_ids


class VariantFilterEvaluator(object):
    """Evaluator for a single scoped expression, e.g. of the form:
        '(position > 5) in ALL(sample1, sample2)'
    """

    def __init__(self, raw_filter_string, ref_genome, scope=None):
        """Constructor.

        Args:
            raw_filter_string: String representing the raw filter.
            ref_genome: ReferenceGenome these variants are relative to.
            scope: Optional FilterScope object which restricts the results
                to the samples according to the semantic setting of the scope.
        """
        # Validation.
        if scope is not None:
            assert isinstance(scope, FilterScope)

        self.raw_filter_string = raw_filter_string
        self.clean_filter_string = raw_filter_string
        self.ref_genome = ref_genome
        self.all_key_map = get_all_key_map(self.ref_genome)
        
        self.variant_caller_common_map = (
                self.ref_genome.get_variant_caller_common_map())
        self.variant_alternate_map = (
                self.ref_genome.get_variant_alternate_map())
        self.variant_evidence_map = (
                self.ref_genome.get_variant_evidence_map())
        self.scope = scope

        # Generator object that provides symbols in alphabetical order.
        self.symbol_maker = symbol_generator()

        # Catch trivial, no filter case.
        if self.clean_filter_string == '':
            self.sympy_representation = ''
        else:
            self._create_symbolic_representation()


    def get_scope_type(self):
        """Returns the scope type.
        """
        if self.scope:
            return self.scope.scope_type
        return None


    def get_scope_sample_id_set(self):
        """Returns the set of sample ids that the scope applies to.
        """
        assert self.scope is not None, (
                "get_scope_sample_id_set() called on none-scoped " +
                "evaluator instance.")
        return self.scope.sample_id_set


    def _create_symbolic_representation(self):
        """Creates a symbolic representation of the query to enable, among
        other things, manipulation with sympy so we can get to disjunctive
        normal form (DNF).
        """
        # Find all the expressions and replace them with symbols.
        self.symbol_to_expression_map = {}

        symbolified_string = self.clean_filter_string
        for regex in [SAMPLE_SCOPE_REGEX, EXPRESSION_REGEX, SET_REGEX,
                GENE_REGEX]:
            symbolified_string = self._symbolify_string_for_regex(
                    symbolified_string, regex)

        self.sympy_representation = boolalg.to_dnf(symbolified_string)


    def _symbolify_string_for_regex(self, start_string, regex):
        """Iterates through the string tokenized by the regex and replaces
        the matching tokens with symbols.
        """
        tokens = regex.split(start_string)
        symbolified_tokens = []
        for token in tokens:
            if regex.match(token):
                symbol = self.symbol_maker.next()
                self.symbol_to_expression_map[symbol] = token
                symbolified_tokens.append(symbol)
            else:
                symbolified_tokens.append(token)
        symbolified_string = ''.join(symbolified_tokens)
        return symbolified_string


    def get_condition_string_for_symbol(self, symbol):
        """Returns the condition string that the symbol replaced.

        Args:
            symbol: A sympy.core.symbol.Symbol or string.

        Returns:
            A string representing a query condition.
        """
        if isinstance(symbol, str):
            symbol_str = symbol
        else:
            symbol_str = str(symbol)
        return self.symbol_to_expression_map[symbol_str]


    def evaluate(self):
        """Evaluates the conditions in the filter string.

        Returns:
            A FilterEvalResult object.
        """
        return self.evaluate_disjunction(self.sympy_representation)


    def evaluate_disjunction(self, disjunction_clause):
        """Evaluate a disjunction (OR) clause.

        Args:
            disjunction_clause: sympy.logic.boolalg.Or clause, or equivalent.
            scope: FilterScope object.

        Returns:
            A FilterEvalResult object.
        """
        # The disjunction may contain a single condition.
        if not isinstance(disjunction_clause, boolalg.Or):
            return self.evaluate_conjunction(disjunction_clause)
        else:
            result = FilterEvalResult(set(), {})
            for conjunction_clause in disjunction_clause.args:
                result |= self.evaluate_conjunction(conjunction_clause)
            return result


    def evaluate_conjunction(self, conjunction_clause):
        """Evaluate a conjunction condition. That is all symbols are ANDed.

        Args:
            disjunction_clause: sympy.logic.boolalg.And clause, or equivalent.
            scope: FilterScope object.

        Returns:
            A FilterEvalResult object.
        """
        # Iterate through the conditions corresponding to the symbols in the
        # clause and either create Q objects out of them, relative to the
        # Variant model, or save them as key-value models to evaluate in memory
        # after the SQL fetch. Order doesn't matter since this is a conjunction
        # clause.

        filter_eval_results = []
        q_list = []
        remaining_triples = []

        if not conjunction_clause == '':
            # The conjunction may contain a single condition.
            if not conjunction_clause.is_Boolean:
                (filter_eval_results, q_list, remaining_triples) = (
                        self._single_symbol_mux(
                                conjunction_clause, filter_eval_results, q_list,
                                remaining_triples))
            else:
                for symbol in conjunction_clause.args:
                    (filter_eval_results, q_list, remaining_triples) = (
                            self._single_symbol_mux(
                                    symbol, filter_eval_results, q_list,
                                    remaining_triples))

        # Combine any filter results so far. These are probably results
        # of sub-clauses that are scoped expressions.
        partial_result = None
        if len(filter_eval_results) > 0:
            partial_result = filter_eval_results[0]
            for result in filter_eval_results[1:]:
                partial_result &= result

        # TODO: Handle evaluating sets for melted view.

        ### Now handle the Q list.

        # NOTE: Previously, we were applying boolean operators
        # among Q objects, but the ORM mapping fails in cases where you
        # want to handle logic among VariantSet membership.
        if len(q_list) > 0:
            variant_list = self.ref_genome.variant_set.filter(q_list[0])
            for q_obj in q_list[1:]:
                variant_list = variant_list.filter(q_obj)
        else:
            variant_list = self.ref_genome.variant_set.all()

        # Make this into a FilterEvalResult object.
        # For now, we just set all the genomes as passing for the conditions
        # so far.
        variant_id_to_metadata_dict = {}
        for variant in variant_list:
            variant_id_to_metadata_dict[variant.id] = {
                'passing_sample_ids': (
                        self.get_sample_id_set_for_variant(variant)),
            }
        q_part_result = FilterEvalResult(set(variant_list),
                variant_id_to_metadata_dict)

        if partial_result is not None:
            partial_result &= q_part_result
        else:
            partial_result = q_part_result

        return self.apply_non_sql_triples_to_query_set(partial_result,
                remaining_triples)


    def _single_symbol_mux(self, symbol, filter_eval_results, q_list,
            remaining_triples):
        """Helper method for evaluating a single symbol.
        """
        result = self.handle_single_symbol(symbol)
        if isinstance(result, FilterEvalResult):
            filter_eval_results.append(result)
        elif isinstance(result, Q):
            q_list.append(result)
        else:
            remaining_triples.append(result)
        return (filter_eval_results, q_list, remaining_triples)


    def handle_single_symbol(self, symbol):
        """Returns one of:
            * A FilterEvalResult object if the symbol represents a scoped
                filter condition.
            * A Django Q object if the symbol represents a condition that
            can be evaluated against the SQL database.
            * A triple of delim, key, value if the condition must be evaluated
                in-memory.
        """
        condition_string = self.get_condition_string_for_symbol(symbol)

        # First check whether this expression is contained within a scope.
        scope_match = SAMPLE_SCOPE_REGEX_NAMED.match(condition_string)
        if scope_match:
            condition_string = scope_match.group('condition')
            scope_type = scope_match.group('scope_type')
            samples_string = scope_match.group('samples')
            sample_ids = FilterScope.parse_sample_ids(samples_string)
            evaluator = VariantFilterEvaluator(condition_string,
                    self.ref_genome, FilterScope(scope_type, sample_ids))
            return evaluator.evaluate()

        # Next, check if this is a set-related expression.
        set_match = SET_REGEX.match(condition_string)
        if set_match:
            return _get_django_q_object_for_set_restrict(condition_string)

        # Next, check if this is a gene-related expression.
        gene_match = GENE_REGEX.match(condition_string)
        if gene_match:
            return _get_django_q_object_for_gene_restrict(condition_string,
                    self.ref_genome)

        # Finally, if here, then this should be a basic, delimiter-separated
        # expression.
        (delim, key, value) = get_delim_key_value_triple(condition_string,
                self.all_key_map)

        # If the key is supported for SQL queries, return the corresponding
        # Q object.
        for key_map in ALL_SQL_KEY_MAP_LIST:
            if key in key_map:
                return _get_django_q_object_for_triple((delim, key, value))

        # Otherwise just return the triple to be evaluated separately.
        return (delim, key, value)


    def apply_non_sql_triples_to_query_set(self, filter_eval_result,
            remaining_triples):
        """Applies the remaining condition triples to the query set and
        returns the trimmed down list.
        """
        # Parse the given FilterEvalResult object.
        variant_list = filter_eval_result.variant_set
        variant_id_to_metadata_dict = (
                filter_eval_result.variant_id_to_metadata_dict)

        for triple in remaining_triples:
            (delim, key, value) = triple

            passing_variant_list = []
            for variant in variant_list:

                # TODO: Currently we are treating alternate keys ENTIRELY as if
                # they were specific to samples. This fine EXCEPT in cases where
                # there exists an INFO value which none of the samples have. The
                # expected result would be to return the variant but no samples,
                # which will not happen if we do it this way.

                # First make sure this is a valid key.
                if not (key in self.variant_caller_common_map or
                        key in self.variant_alternate_map or 
                        key in self.variant_evidence_map):
                    raise ParseError(key, 'Unrecognized filter key.')

                if key in self.variant_caller_common_map:
                    _assert_delim_for_key(self.variant_caller_common_map,
                            delim, key)
                    all_common_data_obj = (
                            variant.variantcallercommondata_set.all())
                    # TODO: Figure out semantics of having more than one common
                    # data object.
                    for common_data_obj in all_common_data_obj:
                        if (_evaluate_condition_in_triple(
                                common_data_obj.as_dict(),
                                self.variant_caller_common_map,
                                triple)):
                            passing_variant_list.append(variant)
                            # No need to update passing sample ids.
                            break

                else: # (if key is per-sample or per-alternate)
                    samples_passing_for_variant = (
                            self.get_samples_passing_for_evidence_or_alternate(
                                    variant, triple))

                    # Determine whether the passing samples qualify this
                    # variant as passing the filter, accounting for scope if
                    # applicable.
                    if len(samples_passing_for_variant) > 0:
                        # Compare the passing results to the scope.
                        scope_type = self.get_scope_type()
                        if (scope_type is None or
                                self.do_passing_samples_satisfy_scope(
                                        scope_type,
                                        samples_passing_for_variant)):
                            passing_variant_list.append(variant)
                            variant_id_to_metadata_dict[variant.id][
                                    'passing_sample_ids'] &= (
                                            samples_passing_for_variant)
                    else:
                        variant_id_to_metadata_dict[variant.id][
                                'passing_sample_ids'] = set()

            # Since we are inside of a conjunction, we only need to check
            # the remaining variants on the next iteration.
            variant_list = passing_variant_list

        return FilterEvalResult(set(variant_list), variant_id_to_metadata_dict)


    def do_passing_samples_satisfy_scope(self, scope_type,
            samples_passing_for_variant):
        assert scope_type in VALID_FILTER_SCOPES, (
                "Unknown scope %s" % scope_type)
        scope_sample_id_set = self.get_scope_sample_id_set()
        if scope_type == FILTER_SCOPE__ALL:
            # All passing sample ids must be in the
            # scope set.
            intersection = (samples_passing_for_variant &
                    scope_sample_id_set)
            if (intersection == scope_sample_id_set):
                return True
        elif scope_type == FILTER_SCOPE__ANY:
            # At least one passing sample id must be in
            # the scope set.
            if len(samples_passing_for_variant &
                    scope_sample_id_set) > 0:
                return True
        elif scope_type == FILTER_SCOPE__ONLY:
            # The passing sample id set must be exactly
            # the scope set.
            if (samples_passing_for_variant ==
                    scope_sample_id_set):
                return True

    def get_samples_passing_for_evidence_or_alternate(self, variant, triple):

        (delim, key, value) = triple
        samples_passing_for_variant = set()

        # Use the appropriate type map if this is a per-alt key.
        if key in self.variant_alternate_map:
            type_map = self.variant_alternate_map
        elif key in self.variant_evidence_map:
            type_map = self.variant_evidence_map
        else:
            raise InputError('Key passed is not in evidence or alternate map.')

        _assert_delim_for_key(type_map, delim, key)

        # Check each VariantEvidence object.
        all_variant_evidence_obj_list = (
                VariantEvidence.objects.filter(
                        variant_caller_common_data__variant=variant))

        for variant_evidence_obj in all_variant_evidence_obj_list:
            data_dict = variant_evidence_obj.as_dict()

            if not data_dict['called']:
                continue

            # For VariantAlternate (per-alt) keys, map the sample's alleles
            # onto items in the list. For instance, if a sample has
            # a genotype of 1/1, then it's items will be the first
            # allele in the -1 list of INFO_EFF_* fields.
            if type_map == self.variant_alternate_map:
                alts = VariantAlternate.objects.filter(
                    variant=variant,
                    variantevidence=variant_evidence_obj)

                data_dict = defaultdict(list)
                for alt in alts:
                    alt_dict = alt.as_dict()
                    [data_dict[k].append(alt_dict[k]) for k in alt_dict.keys()]

            if (_evaluate_condition_in_triple(data_dict, type_map, triple)):
                samples_passing_for_variant.add(
                        variant_evidence_obj.experiment_sample.id)

        return samples_passing_for_variant


    def get_sample_id_set_for_variant(self, variant):
        """Returns the set of all ExperimentSamples ids for which there exists a
        relationship to the Variant
        """
        return set([ve.experiment_sample.id for ve in
                VariantEvidence.objects.filter(
                        variant_caller_common_data__variant=variant)])


###############################################################################
# Helper methods
###############################################################################

def get_per_alt_dict(key, variant, variant_evidence_obj, type_map):
    """Returns a dictionary/type map tuple for one per-alt INFO field.  It is
    relevant to a single variant evidence object, corresponding to all of the
    alternate alleles that the variant evidence object might have.

    For example, given a key of INFO_EFF_SEVERITY, which has the value 
        ['SEVERE','MODERATE']

    corresponding to a SEVERE for alternate allele 1 and a MODERATE severity for
    alternate allele 2.

    Given a variant evidence object where the gt_num is 1/1, the returned
    dictionary would be:

            {'INFO_EFF_SEVERITY':['SEVERE']} 

    If there is an evidene object that has the gt_num  1/2, then the returned
    dictionary would be:

            {'INFO_EFF_SEVERITY':['SEVERE','MODERATE']}

    so that the variant evidence has two chances to match, one if the query
    is looking for SEVERE, and another if it is looking for MODERATE. 

    This function returns a tuple, where the first item is the dictionary
    explained above and the second item is a type map for the particular key
    for used when casting the key outside of its native common_data type map
    dictionary.
    """

    gt_string = variant_evidence_obj.as_dict()['GT']

    # TODO: Could we just ignore phased variants if they are found, and treat
    # them as unphased?
    assert ('|' not in gt_string), (
        'GT string is phased; this is not handled and should never happen...')

    gts = variant_evidence_obj.as_dict()['GT'].split('/')

    # The -1 is because the reference allele is 0, we want the first alt to be 0
    gt_set = set()
    for gt in gts:
        if int(gt) > 0:
            gt_set.add(int(gt)-1)

    key_dict = {key: []}
    data_dict = variant_evidence_obj.variant_caller_common_data.as_dict()

    if key in data_dict:
        evaled_list = data_dict[key]
        for gt in gt_set:
            key_dict[key].append(evaled_list[gt])

    #key_dict[key] = repr(key_dict[key])
    return(key_dict, type_map[key])


def _get_django_q_object_for_triple(delim_key_value_triple):
    """Returns a Django Q object for querying against the SNP model for
    the given key_string.

    Args:
        delim_key_value_triple: A tuple representing a single condition.

    Returns a django Q object.
    """
    assert len(delim_key_value_triple) == 3
    (delim, key, value) = delim_key_value_triple

    # Special handling for != delim.
    if delim == '!=':
        postfix = ''
        maybe_not_prefix = '~'
    else:
        postfix = DELIM_TO_Q_POSTFIX[delim]
        maybe_not_prefix = ''

    eval_string = (maybe_not_prefix + 'Q(' + key + postfix + '=' + '"' + value +
            '"' + ')')
    return eval(eval_string)


def _get_django_q_object_for_set_restrict(set_restrict_string):
    """Returns the Q object for the set restrict.
    """
    match = SET_REGEX_NAMED.match(set_restrict_string)
    variant_set_uid = match.group('sets')
    assert len(variant_set_uid) > 0, (
            "No actual set provided in set filter.")

    q_obj = Q(varianttovariantset__variant_set__uid=variant_set_uid)

    # Maybe negate the result.
    if match.group('maybe_not'):
        return ~q_obj
    return q_obj


def _get_django_q_object_for_gene_restrict(gene_restrict_string, ref_genome):
    """Returns the Q object to limit results to the gene.
    """
    match = GENE_REGEX_NAMED.match(gene_restrict_string)
    gene_label = match.group('gene')

    # Restrict to variants whose position fall between start and end of
    # gene.

    # First attempt, look up the Gene and get positions.
    gene_region = Region.objects.get(
            type=Region.TYPE.GENE,
            reference_genome=ref_genome,
            label=gene_label)

    # Assume gene only has one interval.
    gene_interval = gene_region.regioninterval_set.all()[0]

    # Return a Q object bounding Variants by this position.
    q_obj = (Q(position__gte=gene_interval.start) &
            Q(position__lt=gene_interval.end))
    return q_obj


def _evaluate_condition_in_triple(data_map, type_map, triple, idx=None):
    """Evaluates a condition.

    idx arg loops through all possible items by list index if the data type is a
    list of values (i.e. in the per alternate case - if the data_type map 'spec'
    field is '-1', corresponding to a Number='A' in the vcf). If it is empty,
    then evaluate the data type as one value.

    Idx field calls are  recursive calls from within the function to iterate
    through the list. If any values are true, then the condition returns true.
    """

    (delim, key, value) = triple

    # If this is an INFO field (common_data) and it is a per-alternate field
    # (Number = 'A' in vcf, 'num' == -1 in pyvcf), then match if any of the
    # values is correct. This recursively calls the function with the extra idx
    # field.
    if idx is None and 'num' in type_map[key] and type_map[key]['num'] == -1:
        evaluations = []
        for recurse_idx in range(len(data_map[key])):
            evaluations.append(_evaluate_condition_in_triple(
                    data_map,
                    type_map,
                    triple, 
                    idx=recurse_idx))
        return any(evaluations)
    else:
        cast_type_string = type_map[key]['type']
        if cast_type_string == 'Boolean':
            return _evaluate_boolean_condition(data_map, key, value, idx)
        else:
            casted_value = _cast_value_to_type(value, cast_type_string)
            if idx is not None:
                evaled_list = data_map[key]
                return eval('evaled_list[idx] ' + delim + ' casted_value')
            else:
                return eval('data_map[key] ' + delim + ' casted_value')

def _evaluate_boolean_condition(data_dict, key, value, idx=None):
    """Evaluates a boolean condition.
    """
    VALID_BOOLEAN_TRUE_VALUES = ['True', 'true', 'T', 't']
    VALID_BOOLEAN_FALSE_VALUES = ['False', 'false', 'F', 'f']

    #if data_dict[key] is a dictionary
    if idx is not None:
        init_result = eval(data_dict[key])[idx]
    else:
        init_result = data_dict[key]
    if value in VALID_BOOLEAN_TRUE_VALUES:
        return init_result
    elif value in VALID_BOOLEAN_FALSE_VALUES:
        return not init_result
    else:
        raise ParseError(value, 'Invalid boolean value, use True or False')


def _cast_value_to_type(value, cast_type_string):
    """Return the value casted to the type specified by the cast_type_string,
    as defined in the type maps in the ReferenceGenome object's variant_key_map field
    """
    if cast_type_string == 'Integer':
        return int(value)
    elif cast_type_string == 'Float':
        return float(value)
    elif cast_type_string == 'String':
        return str(value)
    else:
        raise Exception("Unsupported type " + cast_type_string)


def _assert_delim_for_key(type_map, delim, key):
    """Asserts that the delimiter can be evaluated for the type comparison
    specified by the key. Raises a ParseError if not.
    """
    data_type = type_map[key]['type']
    if not delim in TYPE_TO_SUPPORTED_OPERATIONS[data_type]:
        raise ParseError(str(key) + str(delim),
                'Invalid delim for type indicated by key.')

###############################################################################
# Main client method.
###############################################################################

def get_variants_that_pass_filter(filter_string, ref_genome):
    """Takes a complete filter string and returns the variants that pass the
    filter.

    We have a hybrid implementation for field values, where some of the field
    values are actually columns in the respective database tables, while others
    are serialized (pickled) as part of a single key-value string. We do
    filtering over the key-value parts in-memory in the Django server, but on
    a set of objects returned from the SQL database after filtering down by
    the supported DB keys.

    Args:
        filter_string: Query string from the user.
        ref_genome: The ReferenceGenome this is limited to. This is a hard-coded
            parameter of sorts, since at the least we always limit comparisons
            among Variant objects to those that share a ReferenceGenome.

    Returns:
        List of Variant model objects.
    """
    evaluator = VariantFilterEvaluator(filter_string, ref_genome)
    return evaluator.evaluate()
