# Generated by Django 2.0.5 on 2018-05-12 01:24

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('projects', '0005_auto_20180512_0118'),
    ]

    operations = [
        migrations.AddField(
            model_name='projectcategory',
            name='weight',
            field=models.IntegerField(default=1, verbose_name='分类权重'),
        ),
    ]
