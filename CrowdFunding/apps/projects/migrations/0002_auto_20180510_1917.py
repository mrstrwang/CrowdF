# Generated by Django 2.0.5 on 2018-05-10 11:17

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('projects', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='company',
            name='is_approval',
            field=models.BooleanField(default=False, verbose_name='审批状态'),
        ),
        migrations.AddField(
            model_name='project',
            name='is_approval',
            field=models.BooleanField(default=False, verbose_name='审批状态'),
        ),
        migrations.AddField(
            model_name='projectcategory',
            name='is_approval',
            field=models.BooleanField(default=False, verbose_name='审批状态'),
        ),
        migrations.AddField(
            model_name='projectitem',
            name='is_approval',
            field=models.BooleanField(default=False, verbose_name='审批状态'),
        ),
        migrations.AddField(
            model_name='projecttags',
            name='is_approval',
            field=models.BooleanField(default=False, verbose_name='审批状态'),
        ),
    ]
