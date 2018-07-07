from django.db import models
from django.contrib.auth.models import AbstractUser

# Create your models here.


class UserProfile(AbstractUser):
    GENDER_CHOICES = (
        ('male', '男'),
        ('female', '女'),
    )

    ACCOUNT_CHOICES = (
        ('non_choice', '未认证账户'),
        ('company', '商业公司'),
        ('privately', '个体工商户'),
        ('personal', '个人经营'),
        ('government', '政府及非营利组织'),
    )

    ACCOUNT_GROUP = (
        ('member', '会员'),
        ('admin', '管理员'),
    )

    ACCOUNT_TYPE = (
        ('person', '个人'),
        ('company', '企业'),
    )

    account_status = models.CharField(max_length=20, choices=ACCOUNT_CHOICES, default='non_choice', verbose_name='账户类型')
    gender = models.CharField(max_length=6, choices=GENDER_CHOICES, default='male', verbose_name='性别')
    mobile = models.CharField(max_length=11, default='', verbose_name='手机号')
    image = models.ImageField(upload_to='image/%Y/%m', default='img/services-box1.jpg', max_length=100, verbose_name='头像')
    account_group = models.CharField(max_length=20, choices=ACCOUNT_GROUP, default='member', verbose_name='账户组别')
    account_type = models.CharField(max_length=20, choices=ACCOUNT_TYPE, default='person', verbose_name='账户类型')
    id_card_num = models.CharField(max_length=18, default='', verbose_name='身份证')
    id_card_image = models.ImageField(upload_to='users/%Y/%m', default='', verbose_name='身份证照片')
    realname = models.CharField(max_length=20, default='', verbose_name='真实姓名')


    class Meta:
        verbose_name = '用户信息'
        verbose_name_plural = verbose_name

    def __str__(self):
        return self.username


class UserAddress(models.Model):
    user = models.ForeignKey(UserProfile, verbose_name='用户', on_delete='cascade')
    name = models.CharField(max_length=20, verbose_name='收货人')
    mobile = models.CharField(max_length=20, verbose_name='手机')
    address = models.CharField(max_length=100, verbose_name='地址')
    is_default = models.BooleanField(default=False, verbose_name='默认地址')

    class Meta:
        verbose_name = '用户地址'
        verbose_name_plural = verbose_name

    def __str__(self):
        return '{} | {} | {} | {}'.format(self.user, self.name, self.mobile, self.address)