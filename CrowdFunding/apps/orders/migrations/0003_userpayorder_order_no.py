# Generated by Django 2.0.5 on 2018-05-12 03:56

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('orders', '0002_auto_20180510_1910'),
    ]

    operations = [
        migrations.AddField(
            model_name='userpayorder',
            name='order_no',
            field=models.CharField(default='', max_length=100, verbose_name='订单编号'),
        ),
    ]
