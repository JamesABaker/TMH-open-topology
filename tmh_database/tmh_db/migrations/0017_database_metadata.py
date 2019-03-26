# Generated by Django 2.1.7 on 2019-03-26 18:02

from django.db import migrations, models
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('tmh_db', '0016_tmh_hydrophobicity_flexibility'),
    ]

    operations = [
        migrations.CreateModel(
            name='Database_Metadata',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('last_run', models.DateTimeField(default=django.utils.timezone.now)),
                ('last_download', models.DateTimeField(default=django.utils.timezone.now)),
            ],
        ),
    ]
