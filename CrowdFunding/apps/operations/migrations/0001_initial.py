# Generated by Django 2.0.5 on 2018-05-09 16:19

from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('projects', '__first__'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='Advertise',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.CharField(max_length=100, verbose_name='标题')),
                ('weight', models.IntegerField(default=1, verbose_name='广告权重')),
                ('image', models.ImageField(upload_to='advertise/%Y/%m', verbose_name='广告图片')),
                ('project', models.ForeignKey(on_delete='cascade', to='projects.Project', verbose_name='项目')),
            ],
            options={
                'verbose_name': '广告信息',
                'verbose_name_plural': '广告信息',
            },
        ),
        migrations.CreateModel(
            name='Banner',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('is_delete', models.BooleanField(default=False, verbose_name='删除标记')),
                ('create_time', models.DateTimeField(auto_now_add=True, verbose_name='创建时间')),
                ('update_time', models.DateTimeField(auto_now=True, verbose_name='更新时间')),
                ('title', models.CharField(max_length=100, verbose_name='标题')),
                ('image', models.ImageField(upload_to='banner/%Y/%m', verbose_name='轮播图')),
                ('url', models.URLField(verbose_name='访问地址')),
                ('index', models.IntegerField(default=100, verbose_name='访问顺序')),
                ('banner_project', models.ForeignKey(on_delete='cascade', to='projects.Project', verbose_name='众筹项目')),
            ],
            options={
                'verbose_name': '轮播图',
                'verbose_name_plural': '轮播图',
            },
        ),
        migrations.CreateModel(
            name='UserInterest',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('is_delete', models.BooleanField(default=False, verbose_name='删除标记')),
                ('create_time', models.DateTimeField(auto_now_add=True, verbose_name='创建时间')),
                ('update_time', models.DateTimeField(auto_now=True, verbose_name='更新时间')),
                ('project', models.ForeignKey(on_delete='CASCADE', to='projects.Project', verbose_name='关注企业')),
                ('user', models.ForeignKey(on_delete='CASCADE', to=settings.AUTH_USER_MODEL, verbose_name='用户')),
            ],
            options={
                'verbose_name': '用户关注',
                'verbose_name_plural': '用户关注',
            },
        ),
        migrations.CreateModel(
            name='UserSupport',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('is_delete', models.BooleanField(default=False, verbose_name='删除标记')),
                ('create_time', models.DateTimeField(auto_now_add=True, verbose_name='创建时间')),
                ('update_time', models.DateTimeField(auto_now=True, verbose_name='更新时间')),
                ('order_no', models.CharField(max_length=40, verbose_name='订单编号')),
                ('support_nums', models.IntegerField(default=1, verbose_name='回报数量')),
                ('item', models.ForeignKey(on_delete='CASCADE', to='projects.ProjectItem', verbose_name='回报名称')),
                ('project', models.ForeignKey(on_delete='CASCADE', to='projects.Project', verbose_name='项目名称')),
                ('user', models.ForeignKey(on_delete='CASCADE', to=settings.AUTH_USER_MODEL, verbose_name='用户')),
            ],
            options={
                'verbose_name': '用户支持',
                'verbose_name_plural': '用户支持',
            },
        ),
    ]
