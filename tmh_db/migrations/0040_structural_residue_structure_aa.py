# Generated by Django 2.2 on 2019-05-08 10:56

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('tmh_db', '0039_auto_20190508_0959'),
    ]

    operations = [
        migrations.AddField(
            model_name='structural_residue',
            name='structure_aa',
            field=models.CharField(default='X', max_length=1),
        ),
    ]