# Generated by Django 3.1.5 on 2021-04-13 11:39

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('tmh_db', '0033_auto_20210413_1139'),
    ]

    operations = [
        migrations.AddField(
            model_name='structural_residue',
            name='memprotmd_headgroups',
            field=models.BooleanField(null=True),
        ),
        migrations.AddField(
            model_name='structural_residue',
            name='memprotmd_tail',
            field=models.BooleanField(null=True),
        ),
    ]