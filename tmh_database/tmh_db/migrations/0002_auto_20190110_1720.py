# Generated by Django 2.1.5 on 2019-01-10 17:20

from django.db import migrations, models
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('tmh_db', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='proteins',
            name='total_tmh_number',
            field=models.IntegerField(default=None),
        ),
        migrations.AddField(
            model_name='tmhs_locations',
            name='tmh_number',
            field=models.IntegerField(default=0),
            preserve_default=False,
        ),
    ]
