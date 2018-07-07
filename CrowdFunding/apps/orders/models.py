from django.db import models
from utils.Base_model import BaseModel


from users.models import UserProfile, UserAddress
from projects.models import ProjectItem

# Create your models here.


class UserPayOrder(models.Model):
    ORDER_STATUS = (
        ('waiting', '待支付'),
        ('paid', '已支付'),
        ('closed', '订单关闭'),
        ('completed', '已完成'),
    )

    user = models.ForeignKey(UserProfile, verbose_name='用户', on_delete='cascade')
    address = models.ForeignKey(UserAddress, verbose_name='收货地址', on_delete='cascade')
    item = models.ForeignKey(ProjectItem, verbose_name='支持回报', on_delete='cascade')
    item_nums = models.IntegerField(default=1, verbose_name='回报数量')
    freight = models.IntegerField(default=0, verbose_name='配送费用')
    coupon = models.DecimalField(default=0, max_digits=10, decimal_places=2, verbose_name='优惠金额')
    total_pay = models.DecimalField(default=0, max_digits=10, decimal_places=2, verbose_name='总金额')
    need_invoice = models.BooleanField(default=False, verbose_name='是否开票')
    info_invoice = models.CharField(max_length=100, default='个人', verbose_name='发票抬头')
    remarks = models.CharField(max_length=200, null=True, blank=True, verbose_name='备注信息')
    order_time = models.DateTimeField(auto_now_add=True, verbose_name='下单时间')
    pay_time = models.DateTimeField(auto_now_add=True, verbose_name='支付时间')
    order_no = models.CharField(default='', max_length=100, verbose_name='订单编号')
    pay_order_no = models.CharField(default='', max_length=100, verbose_name='付款单号')
    order_status = models.CharField(choices=ORDER_STATUS, default='waiting', max_length=20, verbose_name='订单状态')

    class Meta:
        verbose_name = '订单信息'
        verbose_name_plural = verbose_name

    def __str__(self):
        return '{} | {} | {} | {}'.format(self.user, self.item, self.remarks, self.pay_order_no)