from django.http import Http404
from django.http import HttpResponse
from django.http import HttpResponseBadRequest
from django.http import HttpResponseRedirect
from django.shortcuts import get_object_or_404
from django.shortcuts import render
from django.views.decorators.http import require_POST
from django.conf import settings
from django.views.decorators.csrf import ensure_csrf_cookie
import psycopg2
from boto.s3.connection import Key, S3Connection

import subprocess
import json
import os
import kombu


def get_local_settings_path():
    path = os.path.join(settings.BASE_DIR, "../../genome_designer/conf/local_settings.py")
    assert os.path.isfile(path)
    return path

@ensure_csrf_cookie
def home(request):
    path = get_local_settings_path()
    with open(path) as f:
        code = compile(f.read(), path, 'exec')
    config = {}
    exec code in config

    c = {
        'key': config.get('AWS_SERVER_PUBLIC_KEY', ''),
        'secret': config.get('AWS_SERVER_SECRET_KEY', ''),
        'bucket': config.get('S3_BUCKET', ''),
        'dbname': config.get('DATABASES', {}).get('default', {}).get('NAME', 'genome_designer'),
        'dbname': config.get('DATABASES', {}).get('default', {}).get('NAME', 'genome_designer'),
        'user': config.get('DATABASES', {}).get('default', {}).get('USER', 'genome_designer'),
        'password': config.get('DATABASES', {}).get('default', {}).get('PASSWORD', ''),
        'host': config.get('DATABASES', {}).get('default', {}).get('HOST', ''),
    }
    return render(request, 'base.html', c)

@require_POST
def test_s3(request):
    key = request.POST['key']
    secret = request.POST['secret']
    bucket = request.POST['bucket']

    try:
        S3 = S3Connection(key, secret)
        aws_bucket = S3.get_bucket(bucket)
    except Exception as e:
        print e
        response = {
            'status': 'error',
            'message': "Cannot connect to S3: " + str(e),
        }
        return HttpResponse(json.dumps(response),
            content_type='application/json')

    response = {
        'status': 'success',
    }
    return HttpResponse(json.dumps(response),
        content_type='application/json')


@require_POST
def test_worker(request):
    user = request.POST['user']
    password = request.POST['password']
    host = request.POST['host']
    dbname = request.POST['dbname']

    try:
        connstring = "dbname='%s' user='%s' host='%s' password='%s' port=5432 connect_timeout=5" % (
            dbname, user, host, password)
        conn = psycopg2.connect(connstring)
    except Exception as e:
        print e
        response = {
            'status': 'error',
            'message': "Cannot connect to PostgresSQL: " + str(e),
            'connstring': connstring
        }
        return HttpResponse(json.dumps(response),
            content_type='application/json')

    try:
        connstring = 'amqp://%s:%s@%s:5672//' % (user, password, host)
        conn = kombu.Connection(connstring)
        conn.connect()
    except Exception as e:
        print e
        response = {
            'status': 'error',
            'message': "Cannot connect to RabbitMQ: " + str(e),
            'connstring': connstring
        }
        return HttpResponse(json.dumps(response),
            content_type='application/json')

    response = {
        'status': 'success',
    }
    return HttpResponse(json.dumps(response),
        content_type='application/json')

def save(request):
    user = request.POST['user']
    password = request.POST['password']
    host = request.POST['host']
    dbname = request.POST['dbname']
    key = request.POST['key']
    bucket = request.POST['bucket']
    secret = request.POST['secret']

    db = {
        'default': {
            'ENGINE': 'django.db.backends.postgresql_psycopg2',
            'NAME': dbname,
            'USER': user,
            'PASSWORD': password,
            'HOST': host,
            'PORT': 5432,
        }
    }
    rabbitmq = "amqp://%s:%s@%s:5672//" % (user, password, host)
    config = "DATABASES = %s \n" % repr(db)
    config += "BROKER_URL = %s \n" % repr(rabbitmq)
    config += "S3_BUCKET = %r \n" % repr(bucket)
    config += "AWS_SERVER_PUBLIC_KEY = %r \n" % repr(key)
    config += "AWS_SERVER_SECRET_KEY = %r \n" % repr(secret)
    path = get_local_settings_path()

    with open(path, "w") as f:
        f.write(config)

    os.system("supervisorctl restart celery")

    response = {
        'status': 'success'
    }
    return HttpResponse(json.dumps(response),
        content_type='application/json')

def status(request):
    manage_py_path = os.path.abspath(os.path.join(settings.BASE_DIR, "../../genome_designer/manage.py"))
    commands = ['celery inspect ping',
                'celery status',
                'celery report']

    settings_module = "genome_designer.settings"

    outputs = {}
    for c in commands:
        outputs[c] = os.popen("env DJANGO_SETTINGS_MODULE=%s python %s %s 2>&1" % (
                settings_module, manage_py_path, c)).read()

    response = {
        'status': 'success',
        'outputs': outputs
    }
    return HttpResponse(json.dumps(response),
        content_type='application/json')

