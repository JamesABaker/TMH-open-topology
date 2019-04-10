# Generated by Django 2.2 on 2019-04-10 14:53

from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('tmh_db', '0020_auto_20190410_1408'),
    ]

    operations = [
        migrations.CreateModel(
            name='Funfamstatus',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('submission_key', models.TextField()),
                ('completed_date', models.DateTimeField(default=django.utils.timezone.now)),
                ('protein', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='tmh_db.Protein')),
            ],
        ),
    ]
