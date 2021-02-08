from __future__ import division

import matplotlib
import numpy as np
import pytz
from django.db import models
from django.db.models import Avg
from django.db.models import Case
from django.db.models import Count
from django.db.models import Exists
from django.db.models import F
from django.db.models import Max
from django.db.models import Min
from django.db.models import OuterRef
from django.db.models import Prefetch
from django.db.models import Q
from django.db.models import Subquery
from django.db.models import Sum
from django.db.models import When

from scripts.populate_general_functions import *
from scripts.graphs import *

#def histogram(performance, source, state, x_label, y_label):
tmhs=Tmh_tmsoc.objects.filter(tmh__meta_tmh=True).distinct('pk')


value_list=[]
for i in tmhs:
  value_list.append(i.test_score)

histogram(value_list, "TMSOC", "d", "TMSOC scores", "TMH count")
print(tmhs.count())

#def histogram(performance, source, state, x_label, y_label):
tmhs=Tmh_deltag.objects.filter(tmh__meta_tmh=True).distinct('pk')


value_list=[]
for i in tmhs:
  value_list.append(i.test_score)
print(tmhs.count())

histogram(value_list, "DeltaG", "d", "Delta G scores", "TMH count")
