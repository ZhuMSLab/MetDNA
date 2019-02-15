"""
created by Xin Xiong<onlybelter@outlook.com>, https://github.com/OnlyBelter
"""

import re

from django.conf import settings
from django.db import models
from django.utils.translation import ugettext_lazy as _
from django.contrib.auth.models import (AbstractUser, UserManager)
from django.core.mail import send_mail
from django.core import validators
from django.utils import timezone

# Create your models here.


def user_directory_path(instance, filename):
    # file will be uploaded to
    # MEDIA_ROOT/<username>/<project_hash>/<polarity>/<filename>
    # this is a instance of 'UploadFile' model
    return '{0}/{1}/{2}/{3}'.format(instance.username,
                                    instance.project.project_hash,
                                    instance.polarity,
                                    filename)


class CustomUser(AbstractUser):
    """
    custom user, reference below example
    https://github.com/jonathanchu/django-custom-user-example/blob/master/customuser/accounts/models.py

    # original User class has all I need
    # Just add __str__, not rewrite other field
    - id
    - username
    - password
    - email
    - is_active
    - date_joined
    - method, email_user
    """
    username = models.CharField(_('username'), max_length=30, unique=True,
                                help_text=_('Required. 30 characters or fewer. Letters, numbers and '
                                            '@/./+/-/_ characters'),
                                validators=[
        validators.RegexValidator(re.compile(
            '^[\w.@+-]+$'), _('Enter a valid username.'), 'invalid')
    ])
    full_name = models.CharField(_('full name'), max_length=254, blank=True)
    short_name = models.CharField(_('short name'), max_length=30, blank=True)
    email = models.EmailField(_('email address'), max_length=128, unique=True)
    is_staff = models.BooleanField(_('staff status'), default=False,
                                   help_text=_('Designates whether the user can log into this admin '
                                               'site.'))
    is_active = models.BooleanField(_('active'), default=True,
                                    help_text=_('Designates whether this user should be treated as '
                                                'active. Unselect this instead of deleting accounts.'))
    date_joined = models.DateTimeField(_('date joined'), default=timezone.now)

    objects = UserManager()

    USERNAME_FIELD = 'username'
    REQUIRED_FIELDS = ['email']

    def __str__(self):
        return self.username

    def email_user(self, subject, message, from_email=None):
        """
        Sends an email to this User.
        """
        send_mail(subject, message, from_email, [self.email])


class Project(models.Model):
    user = models.ForeignKey(settings.AUTH_USER_MODEL,
                             related_name='project_id',  # this name will be a field name in CustomUser
                             verbose_name=_('user'),
                             on_delete=models.CASCADE,)
    project_hash = models.CharField(max_length=100,
                                    unique=True,
                                    default='1234567')  # random string, for generating url
    project_name = models.CharField(max_length=100)
    is_active = models.BooleanField(_('active'), default=True,
                                    help_text=_('Designates whether this project should be treated as '
                                                'active. Unselect this instead of deleting accounts.'))

    def __str__(self):
        return self.project_name


class UploadFile(models.Model):
    # this project is an instance of Project
    project = models.ForeignKey(Project,
                                related_name='files_in_project',  # this name will be a field name in Project
                                verbose_name=_('project'),
                                on_delete=models.CASCADE,
                                default=1)
    # current_project = models.CharField(max_length=30, default='null')
    username = models.CharField(max_length=30, default='null')  # set this value in view function
    file_name = models.CharField(max_length=100, blank=True, default='')
    file_type = models.CharField(max_length=20, blank=True)  # ms1_data, ms2_data or sample_info
    stored_url = models.FileField(upload_to=user_directory_path, null=True, blank=True)
    file_size = models.FloatField(max_length=10, blank=True)
    validity = models.CharField(max_length=30, default='1')
    file_info = models.TextField(default='{}')  # default is a serialized empty dict
    polarity = models.CharField(max_length=12, default='NEG')  # POS/NEG/NOT_DEFINED, positive/negative/not defined

    class Meta:
        verbose_name = _('UploadFile')
        verbose_name_plural = _('UploadFiles')
        # https://docs.djangoproject.com/en/dev/ref/models/options/#unique-together
        unique_together = (('project', 'file_name', 'polarity'),)

    def __str__(self):
        return "{0}".format(self.file_name)


class ProjectQueue(models.Model):
    user_id = models.IntegerField(default=1)
    # this project is an instance of Project
    project = models.ForeignKey(Project,
                                # related_name='queue',  # this name will be a field name in Project
                                verbose_name=_('project'),
                                on_delete=models.CASCADE,
                                default=1)
    paras = models.TextField(null=True, blank=True, default='')  # parameters from user input
    status = models.CharField(max_length=10, default='done')
    submit_time = models.DateTimeField(default=timezone.now)
    start_time = models.DateTimeField(default=timezone.now)
    end_time = models.DateTimeField(default=timezone.now)

    def __str__(self):
        return "This project - {0} is {1}".format(self.project.id, self.status)

    class Meta:
        unique_together = (('user_id', 'project'),)
