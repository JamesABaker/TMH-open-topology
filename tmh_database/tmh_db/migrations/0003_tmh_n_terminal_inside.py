# Generated by Django 2.1.7 on 2019-03-21 15:35

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('tmh_db', '0002_tmh_tmh_total_number'),
    ]

    operations = [
        migrations.AddField(
            model_name='tmh',
            name='n_terminal_inside',
            field=models.BooleanField(default=True),
            preserve_default=False,
        ),
    ]
