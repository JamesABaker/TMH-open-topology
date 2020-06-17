# Generated by Django 3.0.4 on 2020-06-11 15:47

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('tmh_db', '0023_auto_20200521_1754'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='Subcellular_location',
            new_name='SubcellularLocation',
        ),
        migrations.RenameField(
            model_name='funfam_residue',
            old_name='e_value',
            new_name='scorecons',
        ),
        migrations.RemoveField(
            model_name='funfam',
            name='uniprot_protein',
        ),
        migrations.AddField(
            model_name='funfam',
            name='superfamily',
            field=models.TextField(default='None'),
        ),
        migrations.AlterField(
            model_name='funfam',
            name='funfam_id',
            field=models.TextField(unique=True),
        ),
        migrations.AlterField(
            model_name='funfam_residue',
            name='funfam',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='tmh_db.Funfam'),
        ),
    ]
