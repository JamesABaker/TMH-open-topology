# Generated by Django 2.2 on 2019-04-10 16:48

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('tmh_db', '0022_auto_20190410_1542'),
    ]

    operations = [
        migrations.AlterField(
            model_name='funfamstatus',
            name='protein',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='tmh_db.Protein', unique=True),
        ),
    ]
