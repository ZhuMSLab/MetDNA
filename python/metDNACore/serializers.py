from rest_framework import serializers
from metDNACore.models import (CustomUser, UploadFile,
                               Project, ProjectQueue)


class UserSerializer(serializers.HyperlinkedModelSerializer):
    project_id = serializers.HyperlinkedRelatedField(many=True,
                                                     view_name='project-detail',
                                                     read_only=True)

    class Meta:
        model = CustomUser
        fields = ('id', 'url', 'username', 'password', 'email',
                  'is_active', 'project_id',)  # 'uploaded_files',
        write_only_fields = ('password', )
        # read_only_fields = ('id', )

    def create(self, validated_data):
        user = CustomUser.objects.create_user(
            username=validated_data['username'],
            email=validated_data['email'],
        )
        user.set_password(validated_data['password'])
        user.save()
        return user

    def update(self, instance, validated_data):
        user = instance
        user.set_password(validated_data['password'])
        user.save()
        return user


class ProjectSerializer(serializers.HyperlinkedModelSerializer):
    files_in_project = serializers.HyperlinkedRelatedField(many=True,
                                                           view_name='uploadfile-detail',
                                                           read_only=True)
    # queue = serializers.HyperlinkedRelatedField(many=False,
    #                                             view_name='project-queue-detail',
    #                                             read_only=True)

    class Meta:
        model = Project
        fields = ('id', 'url', 'project_name', 'project_hash',
                  'user_id', 'files_in_project', 'is_active')


class UploadFileSerializer(serializers.HyperlinkedModelSerializer):
    # user = serializers.ReadOnlyField(source='owner.username')

    class Meta:
        model = UploadFile
        fields = ('id', 'stored_url', 'validity', 'file_type',
                  'polarity', 'project_id', 'file_size', 'file_info')


class ProjectQueueSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = ProjectQueue
        fields = ('id', 'project_id', 'status',
                  'submit_time', 'start_time', 'end_time')
