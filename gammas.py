import sys
sys.path.append('/opt/homebrew/lib/python3.13/site-packages')
from nuclearanalysistools.findGammas import *

# potassium-40, radium-226, and uranium-238


a = AnalyzeGammas(['40K', '226RA', '238U'])
gammas_background = a.findAllGammas()
# gammas_matched = a.matchByGamma(1430.0)
print(gammas_background)

