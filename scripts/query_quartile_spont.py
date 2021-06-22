from scripts.populate_general_functions import *
from scripts.graphs import *
import scipy.stats as stats
import numpy


deltag_objects = Tmh_deltag.objects.filter(tmh_id__meta_tmh=True).distinct("pk")

spont_scores = []
for i in deltag_objects:
    spont_scores.append(i.test_score)

result = numpy.quantile(spont_scores, [0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95])
print(result)
