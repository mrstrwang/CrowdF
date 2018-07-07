from django.db import models

from utils.Base_model import BaseModel
from users.models import UserProfile

# Create your models here.


class ProjectCategory(models.Model):
    name = models.CharField(max_length=10, verbose_name='分类信息')
    desc = models.CharField(max_length=30, verbose_name='分类简介')
    weight = models.IntegerField(default=1, verbose_name='分类权重')
    is_approval = models.BooleanField(default=False, verbose_name='审批状态')

    class Meta:
        verbose_name = '分类信息'
        verbose_name_plural = verbose_name

    def __str__(self):
        return '名称：{} | 简介：{} | 权重：{}'.format(self.name, self.desc, self.weight)


class ProjectTags(models.Model):
    name = models.CharField(max_length=10, verbose_name='一级标签')
    is_approval = models.BooleanField(default=False, verbose_name='审批状态')

    class Meta:
        verbose_name = '项目标签'
        verbose_name_plural = verbose_name

    def __str__(self):
        return self.name


class Company(BaseModel):
    admin = models.ForeignKey(UserProfile, verbose_name='管理用户', on_delete='CASCADE')
    name = models.CharField(max_length=50, unique=True, verbose_name='企业名称')
    desc = models.CharField(max_length=200, verbose_name='企业简介')
    mobile = models.CharField(max_length=20, verbose_name='联系方式')
    image = models.ImageField(upload_to='company/%Y/%m', default='img/p1.jpg', max_length=100, verbose_name='企业图片')
    is_approval = models.BooleanField(default=False, verbose_name='审批状态')

    class Meta:
        verbose_name = '企业名称'
        verbose_name_plural = verbose_name

    def __str__(self):
        return self.name


class Project(BaseModel):
    STATUS_CHOICES = (
        ('readyto', '即将开始'),
        ('funding', '众筹中'),
        ('success', '众筹成功'),
        ('failed', '众筹失败'),
    )

    is_approval = models.BooleanField(default=False, verbose_name='审批状态')
    admin = models.ForeignKey(UserProfile, verbose_name='管理用户', on_delete='CASCADE')
    company = models.ForeignKey(Company, verbose_name='企业名称', on_delete='CASCADE')
    category = models.ForeignKey(ProjectCategory, verbose_name='项目类别', on_delete='CASCADE')
    tags = models.ForeignKey(ProjectTags, verbose_name='项目标签', on_delete='CASCADE')
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='funding', verbose_name='众筹状态')
    name = models.CharField(max_length=100, unique=True, verbose_name='项目名称')
    desc = models.CharField(max_length=200, verbose_name='项目简介')
    interest_nums = models.IntegerField(default=0, verbose_name='关注人数')
    support_nums = models.IntegerField(default=0, verbose_name='支持人数')
    target_fund = models.IntegerField(default=100, verbose_name='目标金额')
    support_fund = models.IntegerField(default=0, verbose_name='已筹资金')
    deadline = models.IntegerField(default=0, verbose_name='筹资天数')
    left_days = models.IntegerField(default=1, verbose_name='剩余天数')
    due_date = models.DateField(auto_now_add=True, verbose_name='到期日期')
    user_desc = models.CharField(max_length=100, verbose_name='自我介绍')
    user_detail = models.CharField(max_length=300, verbose_name='详细自我介绍')
    mobile = models.CharField(max_length=20, verbose_name='联系电话')
    service_mobile = models.CharField(max_length=20, verbose_name='客服电话')
    image_first = models.ImageField(upload_to='project/%Y/%m', max_length=100, verbose_name='项目头图')
    image_detail = models.ImageField(upload_to='project/%Y/%m', max_length=100, verbose_name='项目详情')

    class Meta:
        verbose_name = '众筹项目'
        verbose_name_plural = verbose_name

    def __str__(self):
        return '{}'.format(self.name)


class ProjectItem(BaseModel):
    CATEGORY_CHOICES = (
        ('real', '实物回报'),
        ('virtual', '虚拟物品回报'),
    )

    INVOICE_CHOICES = (
        ('Y', '可开发票'),
        ('N', '不可开发票'),
    )

    project = models.ForeignKey(Project, verbose_name='项目名称', on_delete='CASCADE')
    item_category = models.CharField(choices=CATEGORY_CHOICES, max_length=20, verbose_name='回报类型')
    item_name = models.CharField(max_length=30, verbose_name='回报名称')
    item_value = models.DecimalField(default=1, max_digits=10, decimal_places=2, verbose_name='支持金额')
    item_freight = models.IntegerField(default=0, verbose_name='配送费用')
    item_send_day = models.IntegerField(default=10, verbose_name='发放时间(天)')
    how_to_get = models.CharField(max_length=200, verbose_name='回报获得方式')
    item_nums = models.IntegerField(default=0, verbose_name='回报数量')
    item_image = models.ImageField(upload_to='project/%Y/%m', max_length=100, verbose_name='说明图片')
    support_nums = models.IntegerField(default=0, verbose_name='支持人数')
    up_to_buy = models.IntegerField(default=1, verbose_name='单笔限购')
    invoice = models.CharField(choices=INVOICE_CHOICES, max_length=10, verbose_name='发票')
    is_approval = models.BooleanField(default=False, verbose_name='审批状态')

    class Meta:
        verbose_name = '回报名称'
        verbose_name_plural = verbose_name

    def __str__(self):
        return self.item_name