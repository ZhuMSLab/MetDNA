"""
created by Xin Xiong<onlybelter@outlook.com>, https://github.com/OnlyBelter
"""

import json
import os

from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist
from django.http import JsonResponse
from rest_framework import permissions
from rest_framework import status
from rest_framework import viewsets
from rest_framework.response import Response

from .check_data_validity import check_data
from .models import (CustomUser, UploadFile, Project, ProjectQueue)
from .serializers import (UploadFileSerializer, UserSerializer,
                          ProjectSerializer, ProjectQueueSerializer)
from .send_email import (BASE_MSG_SUBMIT, send_mail)

CURRENT_FRONTEND_URL = 'http://metdna.zhulab.cn'

# Create your views here.


class UserViewSet(viewsets.ModelViewSet):
    permission_classes = (permissions.AllowAny,)
    serializer_class = UserSerializer

    def get_queryset(self):
        # http.get will come here
        queryset = CustomUser.objects.filter(id=-1)
        if self.request.user.id:
            if self.request.user.is_superuser:
                queryset = CustomUser.objects.all()
            else:
                queryset = CustomUser.objects.filter(id=self.request.user.id)
        elif self.request.query_params:
            email = self.request.query_params.get('email', '')
            queryset = CustomUser.objects.filter(email=email)
        return queryset

    def perform_create(self, serializer):
        # http.post will come here
        # print('here has a post')
        # haha
        user_info = dict(self.request.data)
        # print(user_info, type(user_info))
        # import pdb
        # pdb.set_trace()
        serializer.save(
            username=user_info.get('username'),
            password=user_info.get('password'),
            email=user_info.get('email'),
            is_active=user_info.get('is_active')
        )
        # self.perform_create(serializer)
        headers = self.get_success_headers(serializer.data)
        return Response(serializer.data, status=status.HTTP_201_CREATED, headers=headers)

    def perform_update(self, serializer):
        # http.put will come here
        user_info = dict(self.request.data)

        # user_info
        user_info_new = {key: user_info[key][0] for key in user_info.keys()}
        if user_info_new.get('is_active') == 'true':
            user_info_new['is_active'] = True
        else:
            user_info_new['is_active'] = False
        # import pdb
        # pdb.set_trace()
        serializer.save(
            username=user_info_new.get('username'),
            password=user_info_new.get('password'),
            email=user_info_new.get('email'),
            is_active=user_info_new.get('is_active')
        )


class ProjectViewSet(viewsets.ModelViewSet):
    permission_classes = (permissions.AllowAny,)
    serializer_class = ProjectSerializer

    def get_queryset(self):
        queryset = Project.objects.filter(user_id=-1)
        # import pdb
        # pdb.set_trace()
        if self.request.user.id:
            user = self.request.user
            if user.is_superuser:
                queryset = Project.objects.all()
            elif self.request.query_params:
                project_name = self.request.query_params.get('project_name', '')
                if project_name:
                    queryset = Project.objects.filter(user_id=user.id,
                                                      project_name=project_name)
            else:
                queryset = Project.objects.filter(user_id=self.request.user.id)
        return queryset

    def perform_create(self, serializer):
        # we need to get user from HTTP request directly
        # https://stackoverflow.com/a/38162858/2803344
        user = CustomUser(username='test', password='test',
                          email='test@test.com')
        try:
            user = self.request.user
        except ObjectDoesNotExist:
            pass

        # haha

        # print('here has a post')
        project_info = dict(self.request.data)
        if 'is_active' not in project_info:
            project_info['is_active'] = True
        # haha
        if project_info:
            _ = project_info
            # serialize web input data so we can store them into db
            serializer.save(project_hash=_.get('project_hash'),
                            project_name=_.get('project_name'),
                            is_active=_.get('is_active'), user=user)
            headers = self.get_success_headers(serializer.data)
            return Response(serializer.data, status=status.HTTP_201_CREATED, headers=headers)
        else:
            message = 'There is no value was submitted.'
            return JsonResponse(status=404,
                                data={'status': 'false', 'message': message})


class UploadFileViewSet(viewsets.ModelViewSet):
    """
    this viewset needs user and project info, when start a post
    """
    permission_classes = (permissions.AllowAny,)
    serializer_class = UploadFileSerializer
    # queryset = UploadFile.objects.all()

    def get_queryset(self):
        queryset = UploadFile.objects.filter(id=-1)
        # import pdb
        # pdb.set_trace()
        if self.request.user.id:
            if self.request.user.is_superuser:
                queryset = UploadFile.objects.all()
            else:
                queryset = UploadFile.objects.filter(username=self.request.user.username)
        return queryset

    def perform_create(self, serializer):
        # we need to get user from HTTP request directly
        # https://stackoverflow.com/a/38162858/2803344
        user = CustomUser(username='test', password='test',
                          email='test@test.com')
        try:
            user = self.request.user
        except ObjectDoesNotExist as exp:
            print('no user', exp)

        file_info = dict(self.request.data)
        # file_hase is a short string defined by ngx-uploader
        file_hash2type = json.loads(file_info.get('fileId2Type', ['{}'])[0])
        # pol means polarity: pos/neg
        polarity = file_info.get('currentPol', [''])[0].upper()
        file_name2hash = json.loads(file_info.get('fileName2Id', ['{}'])[0])
        project_name = file_info.get('projectName', '0')[0]
        # print(file_info)
        # print(file_hash2type)
        # print(file_name2hash)

        project = Project.objects.filter(id=-1)
        try:
            project = Project.objects.filter(user_id=user.id, project_name=project_name)
            assert len(project) == 1, "project number error!"
            project = project[0]
        except ObjectDoesNotExist as exp:
            print('no project or project name duplicated!', exp)
        # haha
        # print('here has a post')
        my_file = file_info.get('file', '')[0]

        if my_file:
            _ = file_info
            file_name = my_file.name
            file_size = float(my_file.size) / 1024  # kb
            file_hash = file_name2hash.get(file_name, '')
            file_type = file_hash2type.get(file_hash, '')
            # polarity = polarity.get(file_hash, '')
            check_result = {}
            path = os.path.join(settings.MEDIA_ROOT, user.username, project.project_hash)
            if file_type == 'ms1_data' or file_type == 'sample_info':
                # here, need to deal with validity
                check_result = check_data(file_type, my_file, save_path=path, pol=polarity)
            # if file_type == 'ms1_data':
            #     sample_info_file = UploadFile.objects.filter(project_id=project.id, file_type='sample_info')
            #     # import pdb
            #     # pdb.set_trace()
            #     if len(sample_info_file) == 1:
            #         sample_info_file = sample_info_file[0]
            #         print(sample_info_file.polarity, polarity)
            #         if sample_info_file.polarity != polarity:
            #             old_file_path = sample_info_file.stored_url.path
            #             new_file_path = os.path.join(path, polarity, sample_info_file.file_name)
            #             if not os.path.exists(new_file_path):
            #                 shutil.copy(old_file_path, new_file_path)
            # serialize web input data so we can store them into db
            serializer.save(file_name=file_name,
                            file_size=round(file_size, 2),
                            file_type=file_type,  # ms1_data/ms2_data/sample_info
                            stored_url=my_file,
                            username=user.username,
                            project=project,
                            validity=json.dumps(check_result.get('error_code', [0])),
                            file_info=json.dumps(check_result.get('file_info'), {}),
                            polarity=polarity)
            # headers = self.get_success_headers(serializer.data)
            # mes = 'this is a test, 20171017'
            # headers = self.get_success_headers(serializer.data)
            # return Response(serializer.data, status=50000, headers=headers)
        else:
            message = 'There is no value was submitted.'
            return JsonResponse(status=404,
                                data={'status': 'false', 'message': message})

    def perform_destroy(self, instance):
        # http.delete will come here
        # import pdb
        # pdb.set_trace()
        if instance:
            stored_url = str(instance.stored_url)
            file_path = os.path.join(settings.MEDIA_ROOT, stored_url)
            if os.path.exists(file_path):
                os.remove(file_path)
            instance.delete()


class ProjectQueueViewSet(viewsets.ModelViewSet):
    permission_classes = (permissions.AllowAny,)
    serializer_class = ProjectQueueSerializer

    def get_queryset(self):
        queryset = ProjectQueue.objects.filter(id=-1)
        # import pdb
        # pdb.set_trace()
        if self.request.user.id:
            if self.request.user.is_superuser:
                queryset = ProjectQueue.objects.all()
            else:
                queryset = ProjectQueue.objects.filter(user_id=self.request.user.id)
        return queryset

    def perform_create(self, serializer):
        # we need to get user from HTTP request directly
        # https://stackoverflow.com/a/38162858/2803344
        user = CustomUser(username='test', password='test',
                          email='test@test.com')
        try:
            user = self.request.user
        except ObjectDoesNotExist:
            pass
        project = Project.objects.filter(id=-1)
        # haha
        # import pdb
        # pdb.set_trace()
        # print('here has a post')
        project_info = dict(self.request.data)
        project_name = project_info['project_name']
        paras = project_info['paras']
        try:
            project = Project.objects.filter(user_id=user.id, project_name=project_name)
            assert len(project) == 1, "project number error!"
            project = project[0]
        except ObjectDoesNotExist as exp:
            print('no project or project name duplicated!', exp)

        if paras:
            # serialize web input data so we can store them into db
            path = os.path.join(settings.MEDIA_ROOT, user.username, project.project_hash)
            paras['path'] = path
            serializer.save(user_id=user.id,
                            project=project,
                            paras=json.dumps(paras),  # parameters come from user's input
                            status='waiting')
            subject = 'Project {} has submitted!'.format(project.project_name)
            base_msg = BASE_MSG_SUBMIT.format(project_name=project.project_name, faq_link=CURRENT_FRONTEND_URL + '/faq')
            print('subject', subject)
            print('base_msg', base_msg)
            send_mail(to_email=[user.email], subject=subject, message=base_msg)
            # headers = self.get_success_headers(serializer.data)
            # return Response(serializer.data, status=status.HTTP_201_CREATED, headers=headers)
        else:
            message = 'There is no value was submitted.'
            return JsonResponse(status=404,
                                data={'status': 'false', 'message': message})

