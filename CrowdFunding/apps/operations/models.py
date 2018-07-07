from django.db import models
from utils.Base_model import BaseModel

from users.models import UserProfile
from projects.models import ProjectItem, Project

# Create your models here.


class UserSupport(BaseModel):
    user = models.ForeignKey(UserProfile, verbose_name='用户', on_delete='CASCADE')
    project = models.ForeignKey(Project, verbose_name='项目名称', on_delete='CASCADE')
    item = models.ForeignKey(ProjectItem, verbose_name='回报名称', on_delete='CASCADE')
    order_no = models.CharField(max_length=40, verbose_name='订单编号')
    support_nums = models.IntegerField(default=1, verbose_name='回报数量')

    class Meta:
        verbose_name = '用户支持'
        verbose_name_plural = verbose_name

    def __str__(self):
        return '用户：{}， 项目：{}， 订单编号：{}'.format(self.user, self.project, self.order_no)


class UserInterest(BaseModel):
    user = models.ForeignKey(UserProfile, verbose_name='用户', on_delete='CASCADE')
    project = models.ForeignKey(Project, verbose_name='关注企业', on_delete='CASCADE')

    class Meta:
        verbose_name = '用户关注'
        verbose_name_plural = verbose_name

    def __str__(self):
        return '用户：{}， 关注企业：{}'.format(self.user, self.project)


class Banner(BaseModel):
    banner_project = models.ForeignKey(Project, verbose_name='众筹项目', on_delete='cascade')
    title = models.CharField(max_length=100, verbose_name='标题')
    image = models.ImageField(upload_to='banner/%Y/%m', max_length=100, verbose_name='轮播图')
    url = models.URLField(max_length=200, verbose_name='访问地址')
    index = models.IntegerField(default=100, verbose_name='访问顺序')
    is_approval = models.BooleanField(default=False, verbose_name='审批状态')

    class Meta:
        verbose_name = '轮播图'
        verbose_name_plural = verbose_name

    def __str__(self):
        return '{}，当前次序{}'.format(self.title, self.index)


class Advertise(BaseModel):
    title = models.CharField(max_length=100, verbose_name='标题')
    project = models.ForeignKey(Project, verbose_name='项目', on_delete='cascade')
    weight = models.IntegerField(default=1, verbose_name='广告权重')
    image = models.ImageField(upload_to='advertise/%Y/%m', verbose_name='广告图片')
    is_approval = models.BooleanField(default=False, verbose_name='审批状态')

    class Meta:
        verbose_name = '广告信息'
        verbose_name_plural = verbose_name

    def __str__(self):
        return '{} | {} | {}'.format(self.title, self.project, self.weight)


class EmailVerifyCode(models.Model):
    CODE_TYPE = (
        ('register', '注册'),
        ('forget', '忘记密码'),
        ('verify_id', '实名认证'),
    )
    email = models.EmailField(verbose_name='邮箱地址')
    send_time = models.DateTimeField(auto_now_add=True, verbose_name='发送时间')
    verify_code = models.CharField(max_length=100, default=0, verbose_name='邮箱验证码')
    send_type = models.CharField(choices=CODE_TYPE, max_length=20, default='register', verbose_name='验证码类型')

    class Meta:
        verbose_name = '邮箱验证码'
        verbose_name_plural = verbose_name

    def __str__(self):
        return '{} | {} | {} | {}'.format(self.email, self.send_time, self.verify_code, self.send_type)
