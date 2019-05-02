# Generated by Django 2.2 on 2019-05-02 15:59

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('tmh_db', '0034_structural_residue_structure'),
    ]

    operations = [
        migrations.CreateModel(
            name='Binding_residue',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('comment', models.TextField()),
                ('residue', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='tmh_db.Residue')),
            ],
        ),
    ]
