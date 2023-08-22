# Author: Joseph S. Wirth

import copy, logging, os, random, re, time, unittest
from PHANTASM.taxonomy.taxonomyConstruction import constructTaxonomy
from PHANTASM.taxonomy.processLpsn import makeLpsnD
from PHANTASM.taxonomy.Taxonomy import Taxonomy
from PHANTASM.taxonomy.TaxRank import TaxRank
from param import CSV_1, CSV_2, CSV_3, CSV_4
from Bio import Entrez

class TestTaxonomy(unittest.TestCase):
    # initialize the log file
    logging.basicConfig(level=logging.DEBUG, filename="unit_tests.log")
    
    # set entrez.email (required for building Taxonomy objects)
    Entrez.email = 'jwirth@hmc.edu'
    
    # define the temporary filename (used in loadSaveTestSuite)
    TEMP_FN = "temporytestfile.tax"

    # define the max number of children to test (used in runner)
    MAX_NUM_CHILD = 10

    # create lpsnD
    LPSN_D:dict = makeLpsnD(CSV_1, CSV_2, CSV_3, CSV_4)

    # define taxids identified by 16S blastn for each test case
    BLAST_TAX_1  = {1280, 1281, 1282, 1283, 33028, 1285, 1286, 911238, 425474, 1284, 1290, 1276936, 1292, 1293, 1294, 1295, 643214, 1296, 45972, 522262, 569857, 1891096, 1891097, 224030, 246432, 985762, 150056, 46126, 46127, 342451, 94138, 1611836, 29378, 29379, 29380, 1955013, 319942, 71237, 29384, 29385, 214473, 29382, 29388, 155085, 69966, 74703, 1855823, 69968, 198482, 1903056, 198484, 69969, 283734, 1288, 69967, 170573, 74706, 198483, 53344, 201828, 201829, 42858, 1189613, 586733, 70255, 2044912, 145393, 70258, 1141106, 1902580, 28035}
    BLAST_TAX_2  = {859143, 859144, 1504265, 79882, 79883, 220684, 220685, 293386, 79885, 1243664, 293389, 1348623, 1742358, 992288, 2492960, 228899, 673318, 1550385, 119858, 1501239, 1622072, 1808955, 1069116, 762947, 1510468, 580165, 1516104, 483913, 226900, 1707093, 2613334, 1637975, 178774, 1358420, 904295, 64104, 1463404, 2663026, 1349754, 492670, 1926278, 86664, 86668, 33932, 157838, 46224, 324768, 1178786, 283816, 72360, 72361, 218284, 718002, 155322, 2026186, 2026187, 152268, 1150157, 2026190, 2026189, 2026192, 2026193, 2026191, 2026194, 2026188, 2683610, 545502, 1535204, 324854, 137993, 115979, 1285900, 279826, 1208599, 1648923, 285983, 135461, 519977, 356658, 549687, 1347902, 450367, 2478915, 1441095, 1221450, 371036, 363869, 363870, 679261, 363872, 1737571, 531816, 1390, 1037680, 1392, 1396, 1397, 1326968, 1402, 1404, 1405, 1602942, 1408, 421767, 1423, 1428, 1925020, 1925021, 284580, 284581, 574375, 574376, 441769, 1793963, 1452, 1178541, 1178540, 1458, 1467, 189381, 1478, 476102, 189382, 260554, 77777, 889306, 2291674, 866780, 38875, 432607, 302048, 1174504, 209389, 1670641, 96241, 2039284, 264697, 1890302}
    BLAST_TAX_3  = {1307136, 1208324, 644107, 394264, 569882, 639004, 1958944, 647720, 1663016, 1407019, 1917485, 1287727, 2029104, 218673, 1037362, 420403, 1580596, 2546227, 1510456, 1245752, 282683, 393278, 1379903, 1206336, 1775170, 1446476, 1173585, 2494550, 517719, 1396826, 1429085, 2562655, 89184, 1470563, 1679460, 2056810, 74349, 1173614, 999547, 999548, 999549, 999550, 999552, 1155722, 225422, 1246864, 1200281, 1210524, 1200284, 81569, 481446, 1820329, 1499314, 2044597, 1849019, 407234, 2547397, 483013, 1121479, 2578117, 506591, 871651, 871652, 1250539, 396013, 321267, 2499833, 53501, 1695997, 1317121, 580872, 1652492, 2583821, 2583820, 441103, 391948, 1387282, 83219, 996115, 1227549, 571166, 337701, 1844006, 1423144, 2211117, 686389, 321339, 490829, 1280847, 1342301, 1342302, 745311, 2026338, 476529, 1765746, 1156985, 1300350, 1685378, 1685379, 981384, 311180, 314265, 481181, 571298, 505251, 505252, 1492905, 373675, 1123756, 2820523, 1217970, 1010611, 1662395, 393663, 254406, 475083, 42444, 379347, 1294297, 60890, 1294299, 985054, 1715693, 390641, 996342, 555512, 870908}
    BLAST_TAX_4  = {2208, 2209, 523844, 1036677, 2214, 2215, 1036679, 1434110, 38027, 170861, 84592, 109046, 418010, 1715806}
    BLAST_TAX_5  = {523848, 207809, 1702097, 2287}
    BLAST_TAX_6  = {32023, 32024, 983328, 1660067, 824, 1874362, 195, 196, 198, 199, 200, 1031753, 202, 203, 204, 206, 91353, 1244528, 1244531}
    BLAST_TAX_7  = {54291, 1354260, 929813, 1354262, 1563158, 255519, 545, 546, 82979, 548, 82981, 980518, 82983, 150053, 82985, 553, 644651, 556, 480813, 82990, 552, 82989, 82993, 562, 82995, 563, 565, 566, 551989, 82996, 569, 564, 571, 497725, 574, 575, 78398, 573, 1905730, 208962, 580, 81476, 182337, 546367, 1487935, 577, 582, 1175631, 367190, 51288, 53336, 221276, 1177180, 1134687, 1005665, 69218, 1005667, 2529380, 69220, 614, 615, 69222, 158822, 69219, 28151, 881260, 621, 622, 158823, 624, 69224, 623, 158836, 472693, 472694, 472695, 631, 633, 634, 632, 472705, 988801, 158852, 1615494, 1560201, 685706, 263819, 1384589, 244366, 574096, 910996, 300181, 413497, 1296536, 2884250, 1508507, 570012, 158877, 911008, 206499, 65700, 539813, 1510569, 1510570, 413501, 1510573, 1510574, 911022, 160434, 1458355, 642227, 399969, 1161919, 535744, 419007, 1639108, 1300165, 435910, 83655, 2042057, 1341643, 1388748, 61645, 61646, 1972431, 61648, 550, 549, 180435, 61652, 61651, 1028307, 796890, 667127, 1398493, 66270, 66269, 180957, 1926881, 554, 555, 2047724, 68334, 67824, 67825, 67826, 67827, 67828, 67829, 41202, 299767, 2741499, 82988, 1691903, 1667327, 90370, 90371, 1646340, 204038, 204039, 204040, 55208, 204042, 82991, 82992, 566551, 687899, 29471, 1264931, 33060, 2036006, 717607, 1529639, 1737515, 1367852, 29484, 29486, 29483, 47917, 1505588, 58169, 665913, 413499, 1158459, 665914, 413500, 413503, 1324864, 59201, 59202, 59203, 59204, 59205, 1076549, 1076551, 28152, 133448, 92490, 1499973, 137545, 59207, 1076550, 2590031, 1082704, 1073999, 1240404, 138074, 585054, 208224, 371042, 630626, 57706, 1778540, 395631, 357233, 1619313, 35703, 357240, 1509241, 1463164, 935293, 1914243, 2562439, 1201031, 73098, 470931, 470932, 470933, 470934, 39831, 465817, 1477018, 42906, 1109412, 59814, 55207, 2058152, 2675110, 1005994, 55211, 55212, 1336237, 1005996, 1005995, 1006000, 55209, 1006003, 1006004, 1917880, 419257, 451513, 1159613, 1259973, 54736, 368603, 1505757, 693216, 1507808, 1748967, 71656, 28141, 649197, 379893, 1544694, 1544695, 1544696, 1513468, 1449981}
    BLAST_TAX_8  = {316067, 29543, 115783, 2509454, 398767, 2603857, 2603858, 351605, 404380, 1114878, 44671}
    BLAST_TAX_9  = {81409, 53442, 1736, 93718, 1528}
    BLAST_TAX_10 = {186497, 227597, 1151117, 1712654, 227598, 172049, 582419, 69014, 90909, 1128108, 529709, 155321, 163003, 71997, 54077, 71998, 53952, 53953, 195522, 391623, 139207, 523849, 523850, 46539, 272844, 523851, 46540, 110163, 110164, 2261, 2262, 1609559, 593117, 72803, 277988, 187879, 187880, 49899, 49900, 1293037, 71280, 74610, 84597, 54262, 1343739}
    BLAST_TAX_11 = {292563}
    BLAST_TAX_12 = {39688, 1641990, 1641992, 553481, 67081, 758802, 37916, 1520670, 557599, 154654, 1286180, 141349, 2048550, 1286181, 2048551, 228909, 525368, 470074, 470076, 451644, 244292, 152142, 459858, 222805, 487514, 146017, 487521, 146018, 146020, 146021, 33894, 33895, 1273442, 1649257, 56425, 292462, 701041, 701042, 701043, 1041522, 120959, 1209984, 53376, 29314, 53378, 29313, 1962118, 2035343, 1389713, 75922, 115862, 1214102, 1265311, 1265312, 1069220, 431269, 1566886, 144550, 1069221, 444597, 1578165, 85693, 1770, 1341646, 1431246, 1138383, 126673, 912594, 53462, 319705, 319706, 319707, 319709, 187102, 47839, 1775, 1764, 1767, 1648871, 49897, 1858794, 1771, 1772, 1773, 1774, 1768, 1136880, 1777, 1778, 1776, 1780, 1781, 1782, 1783, 1784, 1779, 1786, 1787, 1788, 1789, 1790, 1791, 1792, 1793, 1794, 220927, 1796, 1797, 1798, 39687, 1800, 1799, 1801, 39691, 1804, 39693, 1806, 39695, 1807, 46351, 1810, 1202450, 1809, 39692, 39694, 28445, 212765, 169765, 280871, 43304, 486698, 185642, 722731, 590652, 83262, 117567, 1108812, 2034511, 43348, 1401690, 547163, 370524, 370526, 386911, 216929, 386913, 59750, 398694, 350058, 388459, 1802, 1487726, 561007, 659824, 56689, 153971, 243061, 1260918, 1552759, 83332, 948102, 228230, 1534348, 1534349, 28045, 1203599, 28047, 512402, 59813, 723879, 110505, 560555, 133549, 81858, 354243, 134601, 36809, 36811, 182220, 761804, 258505, 36814, 85968, 404941, 2038731, 36813, 1927124, 36812, 577492, 1795032, 324058, 627680, 84962, 258533, 564198, 44010, 194542, 342002, 188915, 425468, 240125}
    BLAST_TAX_13 = {243232, 579137, 83171, 647171, 2190, 213231, 67760, 1301915, 1069083}
    BLAST_TAX_14 = {2097, 2104, 722438}
    BLAST_TAX_15 = {1166337, 2676234, 172044, 1706000, 413712, 1843237, 1646122, 237611, 237612, 76848, 1541170, 28212, 28213, 2071607, 1664572, 634430, 1355329, 1886787, 53317, 1176646, 1762886, 429134, 2303576, 185949, 185950, 185951, 2116704, 301154, 940133, 152682, 653931, 2707068, 310400, 1141888, 1141890, 1460355, 1141889, 1500313, 1259676, 1779358, 1241761, 330920, 1634989, 511662, 86193, 1219257, 1915074, 298708, 453846, 1462999, 392409, 114404, 40681, 40682, 40683, 1861874, 335605, 1395958, 1395960, 1609977, 99597, 335630, 370959, 1560345, 33051, 563996, 365341, 48935, 344880, 658225, 386874, 363835, 1045317, 1268553, 869719, 463200, 424800, 535906, 535907, 1739108, 1331060, 13688, 13689, 1979270, 93063, 418184, 93064, 977800, 266127, 930197, 1975704, 59803, 1069987, 1344948, 862134, 1344951, 862135, 1123267, 1460676, 397260, 2040274, 68569, 587738, 1813468, 1219043, 1632740, 257003, 1736684, 1763824, 643568, 1736690, 1166323}
    BLAST_TAX_16 = {1017177, 1017181, 104097, 104098, 104099, 104100, 104101, 104102, 108528, 1088869, 110479, 1120919, 112140, 1177712, 1181271, 1203393, 1216886, 1231357, 1231623, 1231624, 1234671, 1234672, 1236500, 1236501, 1236502, 1236503, 1236525, 1260984, 1286186, 1286189, 1286190, 1293412, 1307928, 1307942, 1307948, 1307950, 1307951, 1332824, 1379734, 142834, 146474, 146475, 146476, 1502841, 1502842, 153496, 1552741, 178900, 178901, 182838, 215220, 215221, 220990, 265959, 265960, 272568, 28448, 304077, 318683, 320497, 33995, 33996, 376620, 38307, 38308, 391165, 415421, 431306, 435, 436, 437, 438, 442, 442969, 446692, 481146, 483199, 50715, 634177, 656744, 65958, 65959, 66229, 85325, 864732, 887700, 89584, 91915, 92759, 940265, 940283, 941463}
    BLAST_TAX_17 = {57863, 139, 144, 683292, 1674146, 229155, 373540, 445987, 2761123, 373543, 34095, 521010, 29518, 29519, 100177, 56146, 88916, 664662, 478807, 478174, 314723, 47466, 87162, 42109}
    BLAST_TAX_18 = {1640674, 214856}
    BLAST_TAX_19 = {83554, 83555, 83556, 83557, 83558, 83559, 83560, 1229831, 85991, 813, 331635, 331636, 1405396, 1143323}
    BLAST_TAX_20 = {2052483, 759814, 1848199, 317577, 519440, 1805970, 1299, 57498, 243230, 856736, 1288484, 569127, 2025511, 1768108, 68909, 68910, 745776, 937777, 367282, 1211323, 249403, 249408, 1837379, 1837380, 223556, 1765963, 1452751, 1403343, 1478223, 1867090, 1889238, 1889239, 392408, 693977, 1309407, 1309411, 536443, 1706855, 1182568, 1476583, 55148, 546414, 2136176, 392561, 328692, 394230, 522488, 1007098, 309883, 309884, 1751295}
    BLAST_TAX_21 = {1120986, 49896, 592015, 891968}
    BLAST_TAX_22 = {1460226, 1344003, 334858, 2045452, 1873421, 1873423, 1160719, 103441, 498198, 2072, 75289, 1544730, 675864, 37916, 1046556, 2074, 37919, 75290, 75291, 1544738, 60450, 1286180, 1286181, 236067, 659496, 280618, 2079793, 1678900, 1970232, 39481, 1970234, 451644, 72764, 89154, 1427523, 438851, 2185285, 68170, 152142, 1121362, 286802, 222805, 860246, 81499, 1778267, 118367, 487521, 299618, 464995, 146018, 146021, 72806, 1644129, 292462, 1586287, 1926257, 33907, 33909, 1397367, 1397368, 147577, 33911, 40571, 1280637, 33919, 33920, 1303681, 53378, 1962118, 1522311, 721033, 2662028, 75922, 1214102, 1265311, 283811, 144550, 682667, 1543852, 119981, 280237, 681644, 1670832, 1048753, 1603258, 1530044, 1804989, 64702, 2041025, 479433, 33012, 334542, 341199, 1341646, 47312, 1431246, 630995, 1138383, 1206997, 1750, 53461, 53972, 319705, 319706, 319707, 1752, 1779, 279262, 187102, 29401, 319709, 1755878, 1858794, 1795307, 1772, 629485, 1774, 1095914, 1055466, 1775, 190194, 630515, 228596, 190197, 324851, 84724, 1784, 2016499, 1786, 58107, 2670332, 228602, 1783, 42239, 1792, 2670337, 1794, 54011, 1796, 1788, 2715399, 1802, 39691, 995084, 104205, 1045774, 39693, 86796, 1285901, 1634059, 501011, 501010, 42777, 1818, 415010, 1828, 543527, 43304, 1682732, 99117, 1136941, 103730, 1347891, 103732, 103733, 2069302, 103734, 1843, 103731, 675635, 571190, 590652, 1497396, 83262, 2478913, 2034511, 639313, 102227, 168276, 363865, 43355, 43357, 2020703, 386911, 216929, 386913, 1605990, 57704, 1689449, 1348852, 1716077, 1290095, 182640, 56689, 153971, 1960308, 243061, 1223548, 124285, 400768, 1914755, 2571143, 521096, 239498, 2656650, 714130, 1766802, 1496996, 1735589, 505254, 723879, 59813, 505256, 39846, 1691563, 2654637, 2144174, 76726, 1323732, 913853, 1755582, 359359, 36809, 1387982, 37331, 264148, 556499, 1562584, 2040280, 324058, 2594267, 1670621, 2014, 1730526, 1642469, 433649, 342002, 1365493, 190195, 60920, 398843, 425468, 446462, 1867775}
    BLAST_TAX_23 = {1005928, 1010611, 1031541, 1037362, 1054996, 1123237, 1123247, 1123360, 1123361, 1123755, 1123756, 1125964, 1144477, 1155722, 1156985, 1173584, 1173585, 1200281, 1200284, 1206336, 1217970, 1227549, 1229727, 1247613, 1247867, 1267768, 1267769, 1280831, 1280846, 1280847, 1287727, 1288298, 1294297, 1294299, 1300350, 1312363, 1317121, 1323742, 1342299, 1342301, 1342302, 1356575, 1379903, 1387277, 1391910, 1396826, 1402135, 1407019, 1411902, 1446476, 1449350, 1449351, 1461693, 1461694, 1470562, 1470563, 1481893, 1492771, 1510457, 1510458, 1510460, 1529041, 1537215, 154981, 1563671, 1569283, 1579316, 1580596, 1590614, 1591409, 1608407, 1614693, 1614694, 1632025, 1639690, 1641875, 1652492, 1655543, 1662395, 1679460, 1685378, 1685379, 1695997, 1715691, 1715692, 1715693, 1739787, 1768241, 1775170, 1844006, 1852027, 1852029, 1888910, 1891787, 1914408, 1914409, 1924933, 1940533, 1958944, 1964209, 1977320, 1981892, 2029104, 2072972, 215743, 2170575, 2170577, 2174228, 218672, 218673, 2200892, 2211117, 2211448, 221822, 2250573, 225422, 2259333, 2282382, 2320272, 2340860, 2434, 2483033, 2494550, 2497945, 2499833, 2500159, 252305, 254406, 2545756, 2546227, 2547397, 2547404, 2562279, 2562655, 2576383, 2578117, 2583820, 2583821, 2599291, 2602289, 2605946, 2613965, 266809, 2676438, 2682100, 2691042, 2692236, 2697319, 2716691, 2729614, 2741719, 2747601, 2766981, 2785917, 282683, 2841273, 2865933, 299262, 311180, 314263, 320662, 321267, 321339, 328415, 335975, 337701, 337702, 340021, 344736, 373675, 375761, 379347, 391595, 391948, 393278, 393663, 394264, 396013, 404236, 405746, 420403, 42443, 42444, 441103, 441112, 475083, 481181, 481446, 483013, 488547, 490829, 505252, 505353, 517719, 521013, 53501, 540251, 540747, 55218, 561184, 568098, 569882, 571298, 573024, 576117, 576131, 60137, 60890, 633194, 639004, 644107, 647720, 657014, 670155, 686389, 696760, 696762, 74031, 74348, 74349, 745311, 81569, 83219, 870908, 871651, 871652, 89184, 89187, 928856, 93684, 936889, 981384, 985054, 985055, 996115, 996342, 999547, 999549, 999550, 999552, 999611, 999627, 483324}
    BLAST_TAX_24 = {339658, 645991, 1123369, 592015, 643214, 508460, 187979, 645887, 46469, 49896, 2026194, 248315, 180856, 1397, 1111479, 1096341, 253239, 888060, 437897, 1280706, 180861, 1122184, 926561, 1216062, 551788, 1120987, 33936, 883064, 986075, 2652285, 273068, 866499, 1286, 1502, 1503, 2011011, 572547, 2682456, 2026186, 451755, 220714, 1423, 39778, 661089, 1830138, 471514, 218205, 1408423, 883156, 1428, 2026190, 158847, 281456, 469381, 1296, 29385, 1185412, 885272, 509193, 649639, 509192, 1170224, 1119529, 1220573, 908809, 1868734, 702437, 1925021, 1120971, 1396, 28064, 2655001, 106588, 1120972, 1452, 1516, 580340, 52226, 870242, 2026188, 39777, 1298595, 2026193, 491921, 257985, 1291048, 911092, 246194, 1555112, 2325, 1433, 1123489, 1304874, 2039284, 1053223, 1323375, 553427, 1410665, 520764, 2026189, 584708, 971, 561720, 64104, 574376, 2606906, 679201, 651822, 352165, 447595, 1925020, 1053220, 1392501, 29344, 880478, 1664069, 1005740, 47490, 155322, 1792311, 1405, 1392, 1120975, 178899, 759412, 29350, 927704, 1502943, 1220571, 937334, 296745, 1121468, 36849, 2605789, 1110546, 645512, 891968, 1890302, 1197717, 2754, 81409, 2787712, 1648923, 2064768, 82374, 525903, 189382, 1402, 47055, 500635, 1123249, 1255277, 935948, 180866, 1628085, 1304875, 588581, 926567, 1120997, 2817478, 1528, 103892, 29466, 2026191, 574375, 1053211, 29382, 857293, 930152, 580331, 1712515, 356322, 169759}
    BLAST_TAX_25 = {506, 1499392}

    # here is a sample taxonomy object in table format that can be converted to an actual taxonomy object
    TAX_TABLE = [['69655', 'Sulfurisphaera', 'Sulfurisphaera', 'genus', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '118883', 'False'], ['2100760', 'Saccharolobus', 'Saccharolobus', 'genus', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '118883', 'False'], ['69656', 'Sulfurisphaera ohwakuensis', 'Sulfurisphaera ohwakuensis', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/729/055/GCF_009729055.1_ASM972905v1/GCF_009729055.1_ASM972905v1_genomic.gbff.gz', 'GCF_009729055', 'True', 'TA-1', 'representative genome', '291.0', 'True', '7886111|~|5366651', 'DSM 12421|~|IFO 15161|~|JCM 9065|~|NBRC 15161|~|TA-1|~|DSM12421|~|IFO15161|~|JCM9065|~|NBRC15161', '{}', '69655', 'True'], ['299422', 'Sulfolobus neozealandicus (invalid name)', 'Sulfolobus neozealandicus', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '2284', 'False'], ['43687', 'Metallosphaera sedula', 'Metallosphaera sedula', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'ATCC 51363|~|DSM 5348|~|IFO 15509|~|JCM 9185|~|NBRC 15509|~|TH2|~|ATCC51363|~|DSM5348|~|IFO15509|~|JCM9185|~|NBRC15509', '{}', '41980', 'True'], ['2017961', 'Sulfodiicoccus', 'Sulfodiicoccus', 'genus', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '118883', 'False'], ['73007', 'Sulfolobus yangmingensis', 'Sulfolobus yangmingensis', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'YM1', '{}', '2284', 'False'], ['1938355', 'Sulfolobus mongibelli (invalid name)', 'Sulfolobus mongibelli', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '2284', 'False'], ['518068', 'Acidianus uzoniensis (invalid name)', 'Acidianus uzoniensis', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '12914', 'False'], ['1006005', 'Metallosphaera cuprina', 'Metallosphaera cuprina', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/204/925/GCF_000204925.1_ASM20492v1/GCF_000204925.1_ASM20492v1_genomic.gbff.gz', 'GCF_000204925', 'True', 'Ar-4', 'representative genome', '0.0', 'True', '274058', 'Ar-4|~|CGMCC 1.7082|~|JCM 15769|~|CGMCC1.7082|~|JCM15769', '{}', '41980', 'False'], ['1294262', 'Sulfuracidifex tepidarius', 'Sulfuracidifex tepidarius', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/326/425/GCF_008326425.1_SulIC006_1.0/GCF_008326425.1_SulIC006_1.0_genomic.gbff.gz', 'GCF_008326425', 'True', 'IC-006', 'na', '500.0', 'True', '4539311', 'DSM 104736|~|IC-006|~|JCM 16833|~|DSM104736|~|JCM16833', '{}', '2705406', 'True'], ['1670455', 'Sulfodiicoccus acidiphilus', 'Sulfodiicoccus acidiphilus', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/967/175/GCF_003967175.1_Sacidi_1.0/GCF_003967175.1_Sacidi_1.0_genomic.gbff.gz', 'GCF_003967175', 'True', 'HS-1', 'na', '654.0', 'True', '8119611|~|2181871', 'HS-1|~|InaCC Ar79|~|JCM 31740|~|InaCCAr79|~|JCM31740', '{}', '2017961', 'True'], ['282676', 'Acidianus manzaensis (invalid name)', 'Acidianus manzaensis', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '12914', 'False'], ['1532350', 'Metallosphaera tengchongensis', 'Metallosphaera tengchongensis', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/343/295/GCF_013343295.1_ASM1334329v1/GCF_013343295.1_ASM1334329v1_genomic.gbff.gz', 'GCF_013343295', 'True', 'Ric-A', 'representative genome', '200.0', 'True', '7280851', 'CGMCC 1.12287|~|NBRC 109472|~|Ric-A|~|CGMCC1.12287|~|NBRC109472', '{}', '41980', 'False'], ['207809', 'Sulfolobus tengchongensis (invalid name)', 'Sulfolobus tengchongensis', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '2284', 'False'], ['1111107', 'Metallosphaera yellowstonensis (invalid name)', 'Metallosphaera yellowstonensis', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '41980', 'False'], ['314564', 'Acidianus pozzuoliensis (invalid name)', 'Acidianus pozzuoliensis', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '12914', 'False'], ['47303', 'Sulfuracidifex metallicus', 'Sulfuracidifex metallicus', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/729/515/GCF_009729515.1_ASM972951v1/GCF_009729515.1_ASM972951v1_genomic.gbff.gz', 'GCF_009729515', 'True', 'DSM 6482', 'representative genome', '307.0', 'False', '5376051', 'DSM 6482|~|IFO 15436|~|JCM 9184|~|Kra 23|~|NBRC 15436|~|DSM6482|~|IFO15436|~|JCM9184|~|Kra23|~|NBRC15436', '{}', '2705406', 'False'], ['47304', 'Metallosphaera prunae', 'Metallosphaera prunae', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/005/222/525/GCF_005222525.1_ASM522252v1/GCF_005222525.1_ASM522252v1_genomic.gbff.gz', 'GCF_005222525', 'True', 'Ron 12/II', 'representative genome', '734.0', 'True', '2790751', 'DSM 10039|~|Ron 12/II|~|DSM10039|~|Ron12/II', '{}', '41980', 'False'], ['41673', 'Acidianus brierleyi', 'Acidianus brierleyi', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/201/835/GCF_003201835.2_ASM320183v2/GCF_003201835.2_ASM320183v2_genomic.gbff.gz', 'GCF_003201835', 'True', 'DSM 1651', 'representative genome', '309.0', 'True', '6721991', 'DSM 1651|~|IFO 15269|~|JCM 8954|~|NBRC 15269|~|DSM1651|~|IFO15269|~|JCM8954|~|NBRC15269', '{}', '12914', 'False'], ['41674', 'Stygiolobus', 'Stygiolobus', 'genus', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '118883', 'False'], ['41675', 'Stygiolobus azoricus', 'Stygiolobus azoricus', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/729/035/GCF_009729035.1_ASM972903v1/GCF_009729035.1_ASM972903v1_genomic.gbff.gz', 'GCF_009729035', 'True', 'FC6', 'representative genome', '399.0', 'True', '5366641', 'DSM 6296|~|FC6|~|JCM 9021|~|DSM6296|~|JCM9021', '{}', '41674', 'True'], ['47305', 'Sulfolobus thuringiensis (invalid name)', 'Sulfolobus thuringiensis', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '2284', 'False'], ['43080', 'Sulfolobus islandicus (invalid name)', 'Sulfolobus islandicus', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '2284', 'False'], ['1702097', 'Saccharolobus caldissimus', 'Saccharolobus caldissimus', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/886/315/GCF_020886315.1_ASM2088631v1/GCF_020886315.1_ASM2088631v1_genomic.gbff.gz', 'GCF_020886315', 'True', 'JCM 32116', 'representative genome', '517.0', 'True', '11456901', 'HS-3|~|InaCC Ar80|~|JCM 32116|~|InaCCAr80|~|JCM32116', '{}', '2100760', 'False'], ['111955', 'Sulfurisphaera tokodaii', 'Sulfurisphaera tokodaii', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/205/GCF_000011205.1_ASM1120v1/GCF_000011205.1_ASM1120v1_genomic.gbff.gz', 'GCF_000011205', 'True', '7', 'representative genome', '0.0', 'True', '32368', '7|~|DSM 16993|~|JCM 10545|~|NBRC 100140|~|DSM16993|~|JCM10545|~|NBRC100140', '{}', '69655', 'False'], ['452949', 'Sulfolobus beitou (invalid name)', 'Sulfolobus beitou', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '2284', 'False'], ['452950', 'Sulfolobus vallisabyssus (invalid name)', 'Sulfolobus vallisabyssus', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '2284', 'False'], ['2049879', 'Sulfurisphaera javensis', 'Sulfurisphaera javensis', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'InaCC Ar81|~|JCM 32117|~|KD-1|~|InaCCAr81|~|JCM32117', '{}', '69655', 'False'], ['312539', 'Acidianus sulfidivorans', 'Acidianus sulfidivorans', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/201/765/GCF_003201765.2_ASM320176v2/GCF_003201765.2_ASM320176v2_genomic.gbff.gz', 'GCF_003201765', 'True', 'JP7', 'representative genome', '606.0', 'True', '6721981', 'DSM 18786|~|JCM 13667|~|JP7|~|DSM18786|~|JCM13667', '{}', '12914', 'False'], ['2824673', 'Stygiolobus caldivivus', 'Stygiolobus caldivivus', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/704/315/GCF_019704315.1_ASM1970431v1/GCF_019704315.1_ASM1970431v1_genomic.gbff.gz', 'GCF_019704315', 'True', 'KN-1', 'representative genome', '896.0', 'True', '10752551', 'JCM 34622|~|KCTC 4293|~|KN-1|~|JCM34622|~|KCTC4293', '{}', '41674', 'False'], ['118883', 'Sulfolobaceae', 'Sulfolobaceae', 'family', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '2281', 'False'], ['269667', 'Acidianus convivator (invalid name)', 'Acidianus convivator', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '12914', 'False'], ['146920', 'Acidianus tengchongensis (invalid name)', 'Acidianus tengchongensis', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '12914', 'False'], ['563177', 'Acidianus hospitalis (invalid name)', 'Acidianus hospitalis', 'species', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '12914', 'False'], ['2283', 'Acidianus ambivalens', 'Acidianus ambivalens', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/428/885/GCF_009428885.1_ASM942888v1/GCF_009428885.1_ASM942888v1_genomic.gbff.gz', 'GCF_009428885', 'True', 'DSM 3772', 'na', '30.0', 'False', '5091851', 'CIP 104912|~|DSM 3772|~|JCM 9191|~|Lei 10|~|CIP104912|~|DSM3772|~|JCM9191|~|Lei10', '{}', '12914', 'False'], ['2284', 'Sulfolobus', 'Sulfolobus', 'genus', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '118883', 'True'], ['2285', 'Sulfolobus acidocaldarius', 'Sulfolobus acidocaldarius', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/472/365/GCF_028472365.1_ASM2847236v1/GCF_028472365.1_ASM2847236v1_genomic.gbff.gz', 'GCF_028472365', 'True', 'DSM 639', 'na', '419.109', 'True', '15805891', 'ATCC 33909|~|DSM 639|~|IFO 15157|~|JCM 8929|~|NBRC 15157|~|ATCC33909|~|DSM639|~|IFO15157|~|JCM8929|~|NBRC15157', '{}', '2284', 'True'], ['2286', 'Saccharolobus shibatae', 'Saccharolobus shibatae', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/175/345/GCF_019175345.1_ASM1917534v1/GCF_019175345.1_ASM1917534v1_genomic.gbff.gz', 'GCF_019175345', 'True', 'B12', 'representative genome', '40.0', 'True', '10472031', 'ATCC 51178|~|B 12|~|DSM 5389|~|IFO 15437|~|JCM 8931|~|NBRC 15437|~|ATCC51178|~|B12|~|DSM5389|~|IFO15437|~|JCM8931|~|NBRC15437', '{}', '2100760', 'False'], ['2287', 'Saccharolobus solfataricus', 'Saccharolobus solfataricus', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/079/115/GCF_900079115.1_SSOP1/GCF_900079115.1_SSOP1_genomic.gbff.gz', 'GCF_900079115', 'True', 'P1', 'na', '65.0', 'True', '732271', 'ATCC 35091|~|DSM 1616|~|IFO 15331|~|JCM 8930|~|NBRC 15331|~|P1|~|ATCC35091|~|DSM1616|~|IFO15331|~|JCM8930|~|NBRC15331', '{}', '2100760', 'True'], ['79601', 'Metallosphaera hakonensis', 'Metallosphaera hakonensis', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/201/675/GCF_003201675.2_ASM320167v2/GCF_003201675.2_ASM320167v2_genomic.gbff.gz', 'GCF_003201675', 'True', 'HO1-1', 'representative genome', '647.0', 'True', '6721971', 'ATCC 51241|~|DSM 7519|~|HO1-1|~|IAM 14250|~|JCM 8857|~|ATCC51241|~|DSM7519|~|IAM14250|~|JCM8857', '{}', '41980', 'False'], ['12914', 'Acidianus', 'Acidianus', 'genus', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '118883', 'False'], ['12915', 'Acidianus infernus', 'Acidianus infernus', 'species', 'False', 'True', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/729/545/GCF_009729545.1_ASM972954v1/GCF_009729545.1_ASM972954v1_genomic.gbff.gz', 'GCF_009729545', 'True', 'DSM 3191', 'representative genome', '348.0', 'False', '5376081', 'DSM 3191|~|IFO 15270|~|JCM 8955|~|NBRC 15270|~|So4a|~|DSM3191|~|IFO15270|~|JCM8955|~|NBRC15270', '{}', '12914', 'True'], ['41980', 'Metallosphaera', 'Metallosphaera', 'genus', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '118883', 'False'], ['2705406', 'Sulfuracidifex', 'Sulfuracidifex', 'genus', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', '118883', 'False'], ['2281', 'Sulfolobales', 'Sulfolobales', 'order', 'False', 'True', 'None', 'None', 'False', 'None', 'None', '0.0', 'False', '', 'None', '{}', 'None', 'False']]

    # define functions so that each test's time will be printed
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        duration = time.time() - self.startTime
        duration = int(round(duration))

        MSG_1 = "\n** ran test case in "
        MSG_2 = " seconds **\n"

        print(MSG_1 + str(duration) + MSG_2)
        logging.info(MSG_1[1:] + str(duration) + MSG_2)
        time.sleep(1)

        # remove the temp file if it is present
        if os.path.exists(TestTaxonomy.TEMP_FN):
            os.remove(TestTaxonomy.TEMP_FN)



    # Taxonomy object construction
    def buildTaxO(blasttax:set) -> Taxonomy:
        """ buildTaxO:
                Accepts a set of blast taxids as input. Constructs and returns
                a Taxonomy object.
        """
        # constants for printing
        MSG_1 = "\n** building test case **"
        MSG_2 = "** built test case in "
        MSG_3 = " seconds **\n"
        
        # required for building the Taxonomy objects
        Entrez.email = 'jwirth@cdc.gov'

        # log start time and print message
        start = time.time()
        print(MSG_1)

        # build object
        taxO = constructTaxonomy(blasttax, saveTax=False)

        # log total time, print message, and return
        end = time.time()
        duration = int(round((end-start)))
        print(MSG_2 + str(duration) + MSG_3)
        logging.info(MSG_2[1:] + str(duration) + MSG_3[:-1])
        return taxO




    # test suites (majority of all tests defined here)
    def taxonomyConstructionTestSuite(self, taxO:Taxonomy) -> None:
        """ taxonomyConstructionTestSuite:
                Accepts a Taxonomy object as input. Ensures that the object ma-
                tches the expected format and that it is consistent. Does not
                return.
        """
        # constants
        SPECIES_FORMAT = r"^\S+ \S+$"
        INVALID_STR = ' (invalid name)'
        PLACEHOLD_FIND = r'^\[(\S+); .+$'
        PLACEHOLD_REPL = r'\1'
        INVALID_GREP = r' \(invalid name\)$'
        QUOTE_GREP = r'"'
        HIGHER_TAX_FORMAT = r"^\S+$"
        EMPTY_STR = ''

        # make sure it is a Taxonomy object
        self.assertIsInstance(taxO, Taxonomy)

        # check non-Taxonomy member variables whose types should not vary
        self.assertIsInstance(taxO.taxid, str)
        self.assertIsInstance(int(taxO.taxid), int)
        self.assertIsInstance(taxO.sciName, str)
        self.assertIsInstance(taxO.ncbiName, str)
        self.assertIsInstance(taxO.rank, TaxRank)
        self.assertIsInstance(taxO.descendantsD, dict)
        self.assertIsInstance(taxO.isExternal, bool)
        self.assertIsInstance(taxO.ncbiNameIsCorrect, bool)
        self.assertIsInstance(taxO.assemblyFromType, bool)
        self.assertIsInstance(taxO.assCoverage, float)
        self.assertIsInstance(taxO.assIsComplete, bool)
        self.assertIsInstance(taxO.allAssIds, list)
        self.assertIsInstance(taxO.taxidSynonyms, dict)

        # make sure '0' is not being used as a taxid
        self.assertFalse(taxO.getRoot().getDescendantByTaxId(0))

        # check that all scientific names have the expected format
        root = taxO.getRoot()
        curRank = copy.deepcopy(root.rank)

        # we will exit the loop after checking species
        while True:
            # get the list of descendants to evaluate
            if curRank == root.rank:
                allDesc = [root]
            else:
                allDesc = root.getDescendantsAtRank(curRank, list)

            # for each descendant taxon
            desc:Taxonomy
            for desc in allDesc:
                # unpack any invalid names
                temp = re.sub(INVALID_GREP, EMPTY_STR, desc.sciName)
                temp = re.sub(QUOTE_GREP, EMPTY_STR, temp)

                # unpack any placeholder names
                temp = re.sub(PLACEHOLD_FIND, PLACEHOLD_REPL, temp)

                if curRank == Taxonomy.SPECIES:
                    # species names should be exactly two words
                    self.assertTrue(bool(re.search(SPECIES_FORMAT, temp)))

                    # un-modify any invalid names (remove invalid and quotes)
                    parName = desc.parent.sciName
                    parName = re.sub(QUOTE_GREP, EMPTY_STR, parName)
                    parName = re.sub(INVALID_GREP, EMPTY_STR, parName)

                    # validly named species should contain their parental genus
                    if INVALID_STR not in desc.sciName:
                        self.assertTrue(bool(re.search(parName, temp)))

                else:
                    
                    # non-species should be exactly one word long
                    self.assertTrue(bool(re.search(HIGHER_TAX_FORMAT, temp)))
            
            # keep checking lower ranks until species is reached
            if curRank > Taxonomy.SPECIES: curRank.decrement()
            else: break

        # confirm consistency
        ### checks parents are listed as the parent of their children
        ### checks children are listed as the children of their parent
        ### checks that taxonomic ranks are distributed in a rational way
        ### checks that there are no shared taxids (keys) within the object
        ### checks that the type material is always also a child
        ### checks that scientific names are not shared within the object
        self.assertTrue(taxO._isConsistent())

        # confirm the root's parent is None
        root = taxO.getRoot()
        self.assertEqual(root.parent, None)

    def equalityCopyTestSuite1(self, taxO:Taxonomy) -> None:
        """ equalityCopyTestSuite1:
                Accepts a Taxonomy object as input. Evaluates the expected fun-
                ctionality for Taxonomy.copy() and Taxonomy.__eq__(). Specific-
                ally, it checks that a copied object is equal to the original,
                and that multiple, subsequent copies maintain this feature.
                Does not return.
        """
        # check single copy
        self.assertEqual(taxO, taxO.copy())
        self.assertEqual(taxO.copy(), taxO.copy())

        # check double copy
        self.assertEqual(taxO, taxO.copy().copy())
        self.assertEqual(taxO.copy(), taxO.copy().copy())

    def equalityCopyTestSuite2(self, taxO:Taxonomy) -> None:
        """ equalityCopyTestSuite2:
                Accepts a Taxonomy object as input. Evaluates the expected fun-
                ctionality for Taxonomy.copy() and Taxonomy.__eq__(). Specific-
                ally, it checks that copy is a deepcopy and that modifying the
                copy does not modify the original object or any other unmodifi-
                ed copies. In this test suite, the modified member variable is
                NOT a Taxonomy object. Does not return.
        """
        # check modifying a non-Taxonomy member variable
        copy1 = taxO.copy()  # will be modified
        copy2 = taxO.copy()  # will remain unmodified

        # if the object's rank is greater than species
        if copy1.rank > Taxonomy.SPECIES:
            # set 'taxon' to a random descendant species
            allSpecies = list(copy1.getDescendantsAtRank('species').values())
            taxon:Taxonomy = random.choice(allSpecies)
        
        # otherwise the object is a species; only one choice to make
        else:
            taxon:Taxonomy = copy1
        
        # modify the selected species object
        taxon.assemblyFtp = 'touch'

        # the original is different from a modified copy
        self.assertNotEqual(taxO, copy1)

        # an unmodified copy is different from a modified copy
        self.assertNotEqual(copy2, copy1)

        # the original is identical to an unmodified copy
        self.assertEqual(taxO, copy2)

    def equalityCopyTestSuite3(self, taxO:Taxonomy) -> None:
        """ equalityCopyTestSuite3:
                Accepts a Taxonomy object as input. Evaluates the expected fun-
                ctionality for Taxonomy.copy() and Taxonomy.__eq__(). Specific-
                ally, it checks that copy is a deepcopy and that modifying the
                copy does not modify the original object or any other unmodifi-
                ed copies. In this test suite, the modified member variable is
                the typeMaterial Taxonomy object. Does not return.
        """
        # check modifying a Taxonomy member variable
        copy1 = taxO.copy()  # will be modified
        copy2 = taxO.copy()  # will remain unmodified

        # for ranks higher than species, pick a random descendant from any rank
        if copy1.rank > Taxonomy.SPECIES:
            allDesc = copy1.getAllDescendants(list)
            taxon = random.choice(allDesc)
        
        # if the object is a species, it has no descendants, so just use it
        else:
            taxon = copy1
        
        # modify the typeMaterial field ...
        if taxon.rank > Taxonomy.SPECIES:
            # ... ranks higher than species use Taxonomy objects (or None)
            taxon.typeMaterial = Taxonomy(42, "fakeTax", \
                                          taxon.rank.getRankBelow())
        else:
            # ... species always use a list (can be empty, but never None)
            taxon.typeMaterial = ['fake', 'strain', 'list']

        # the original is different from a modified copy
        self.assertNotEqual(taxO, copy1)

        # an unmodified copy is different from a modified copy
        self.assertNotEqual(copy2, copy1)

        # the original is identical to an unmodified copy
        self.assertEqual(taxO, copy2)
    
    def equalityCopyTestSuite4(self, taxO:Taxonomy) -> None:
        """ equalityCopyTestSuite4:
                Accepts a Taxonomy object as input. Evaluates the expected fun-
                ctionality of Taxonomy.copy() and Taxonomy.__eq__(). Specifica-
                lly, it checks that copy is a deepcopy and that modifying the
                copy does not modify the original object or any other unmodifi-
                ed copies. In this test suite, the modified member variable is
                the parent Taxonomy object. The equality operator specifically
                DOES NOT evaluate the lineage, but DOES evaluate the descendant
                Taxonomy objects. In other words, "the sins of the father" can-
                not be inherited. Does not return.
        """
        # only proceed if their is a parent to modify (skip the root)
        if taxO != taxO.getRoot():
            # verify that modifying a parent DOES NOT make objects inequivalent
            copy1 = taxO.copy()  # will be modified
            copy2 = taxO.copy()  # will remain unmodified
        
            # modify the parental object
            copy1.parent.assCoverage = -42  # impossible value

            # despite modifying copy1's parent, all objects should be equal
            self.assertEqual(taxO, copy1)
            self.assertEqual(copy2, copy1)
            self.assertEqual(taxO, copy2)

    def loadSaveTestSuite(self, taxO:Taxonomy) -> None:
        """ loadSaveTestSuite:
                Accepts a Taxonomy object as input. Evaluates the expected fun-
                ctionality of Taxonomy.save() and Taxonomy.load(). Affirms that
                an object can be saved, that an identical object can be loaded
                from the save file, and that the process did not modify the 
                Taxonomy objects in any way.
        """
        # make a copy as a way to compare to the original state
        copy1 = taxO.copy()

        # save the object, load a new Taxonomy object from file, and then rm
        taxO.save(TestTaxonomy.TEMP_FN)
        load:Taxonomy = Taxonomy.load(TestTaxonomy.TEMP_FN)
        os.remove(TestTaxonomy.TEMP_FN)

        # loaded object (handled by the root) is equal to the original
        self.assertEqual(taxO.getRoot(), load.getRoot())

        # the original object has not been modified during the save/load process.
        self.assertEqual(taxO, copy1)

    def assemblyTestSuite(self, taxO:Taxonomy) -> None:
        """ assemblyTestSuite:
                Accepts a Taxonomy object as input. Evaluates the various memb-
                er variables that store NCBI Assembly data.
        """
        # get a list of all the species
        speciesL = taxO.getDescendantsAtRank('species', list)

        # for each species
        speO:Taxonomy
        for speO in speciesL:
            # ensure that an assembly does not exist if no strain present
            if speO.assemblyStrain is None:
                self.assertIsNone(speO.assemblyFtp)

            # if an assembly exists ...
            if speO.assemblyFtp is not None:
                # ... ensure that a strain is designated
                self.assertIsInstance(speO.assemblyStrain, str)

                # ... ensure that assemblyFromType matches expectations
                if speO.assemblyFromType:
                    self.assertIn(speO.assemblyStrain, speO.typeMaterial)
                else:
                    self.assertNotIn(speO.assemblyStrain, speO.typeMaterial)
                
                # ... ensure that the accession number matches the ftp path
                self.assertIn(speO.assemblyAccn, speO.assemblyFtp)

    def typeMaterialTestSuite(self, taxO:Taxonomy) -> None:
        """ typeMaterialTestSuite:
                Accepts a Taxonomy object as input. Evaluates the typeMaterial
                field for the expected format. Some of the tests for higher ra-
                nks should have overlap with the consistencyRootTestSuite if
                Taxonomy._isConsistent() is working as expected. Does not retu-
                rn.
        """
        # type material for species is a list (empty ok) of strings
        if taxO.rank == Taxonomy.SPECIES:
            self.assertIsInstance(taxO.typeMaterial, list)
            for strain in taxO.typeMaterial:
                self.assertIsInstance(strain, str)
        
        # type material for higher ranks ...
        else:
            # ... is Nonetype or Taxonomy
            self.assertIn(type(taxO.typeMaterial), {type(None), Taxonomy})

            if taxO.typeMaterial is not None:
                # ... is also a child of the calling object
                self.assertIn(taxO.typeMaterial, taxO.getChildren(set))

                # ... lists the calling object as its parent
                self.assertEqual(taxO, taxO.typeMaterial.parent)

                # ... is exactly one rank below the calling object
                self.assertEqual(taxO.rank.getRankBelow(), \
                                 taxO.typeMaterial.rank)
            
            # recurse on all children
            for child in taxO.getChildren(set):
                self.typeMaterialTestSuite(child)

    def lpsnTestSuite1(self, taxO:Taxonomy) -> None:
        """ lpsnTestSuite1:
                Accepts a Taxonomy object as input. Tests if the LPSN data were
                fully amalgamated into the Taxonomy object. Specifically, tests
                that all valid names are in lpsnD and that the typeMaterial and
                parent fields were populated correctly. Does not return.
        """
        # constants
        TYPE = 'type'
        PARENT = 'parent'
        INVALID_STR = "(invalid name)"
        PLACEHOLDER_STR = "; not assigned to "

        # use a more compact name
        lpsnD = TestTaxonomy.LPSN_D

        # get the keys for accessing lpsnD
        rankKey = str(taxO.rank)

        # if the name is not marked as invalid ...
        if INVALID_STR not in taxO.sciName:
            # ... and it is a species ...
            if taxO.rank == Taxonomy.SPECIES:
                # ... then every strain available for import should be present
                for strain in lpsnD[rankKey][taxO.sciName][TYPE]:
                    self.assertIn(strain, taxO.typeMaterial)
            
            # ... and it is a rank higher than species ...
            else:
                # ... then it should match lpsnD at the typeMaterial field
                if taxO.typeMaterial is not None:
                    self.assertEqual(taxO.typeMaterial.sciName, \
                        lpsnD[rankKey][taxO.sciName][TYPE])
                
                # ... then it should match lpsnD at the parent field ...
                if taxO.parent is not None:


                    lpsnParent = lpsnD[rankKey][taxO.sciName][PARENT]
                    # ... assuming that the name is not an LPSN placeholder
                    if PLACEHOLDER_STR not in lpsnParent:
                        self.assertEqual(taxO.parent.sciName, lpsnParent)

                    
                    # ... unless LPSN is using a placeholder ...
                    else:
                        # ... then it should be equal or different correct name
                        parRankKey = str(taxO.parent.rank)
                        self.assertTrue(taxO.parent.sciName == lpsnParent or \
                              taxO.parent.sciName in lpsnD[parRankKey].keys())
        
        # recurse on the children
        for child in taxO.getChildren(set):
            self.lpsnTestSuite1(child)

    def lpsnTestSuite2(self, taxO:Taxonomy) -> None:
        """ lpsnTestSuite2:
                Accepts a Taxonomy object as input. Tests if the LPSN data were
                fully amalgamated into the Taxonomy object. Specifically, tests
                that there are no knonwn synonyms stored at the sciName field 
                anywhere within the object. Does not return.
        """
        # constants
        SYN_KEYS = {'species_syn', 'genus_syn', 'family_syn', 'order_syn'}

        # ensure that no synonyms are present in the taxonomy object
        for synKey in SYN_KEYS:
            for synName in TestTaxonomy.LPSN_D[synKey]:
                self.assertFalse(taxO.getDescendantBySciName(synName))
    
    def lpsnTestSuite3(self, taxO:Taxonomy) -> None:
        """ lpsnTestSuite3:
                Accepts a Taxonomy object as input. Tests if the LPSN data were
                fully amalgamated into the Taxonomy object. Specifically, tests
                that the names marked as invalid are completely absent from the
                LPSN data. Does not return.
        """
        # constants
        INVALID_STR = "(invalid name)"
        SYN_SUFFIX  = "_syn"

        # if a name is marked as invalid
        if INVALID_STR in taxO.sciName:
            # obtain the relevant dictionaries to query
            rankD:dict = TestTaxonomy.LPSN_D[str(taxO.rank)]
            synD:dict  = TestTaxonomy.LPSN_D[str(taxO.rank)+SYN_SUFFIX]
            
            # the ncbiName should be within the sciName
            self.assertIn(taxO.ncbiName, taxO.sciName)

            # the ncbiName should be absent from both lpsn dictionaries
            self.assertNotIn(taxO.ncbiName, rankD.keys())
            self.assertNotIn(taxO.ncbiName, synD.keys())
        
        # recurse on the children
        for child in taxO.getChildren(set):
            self.lpsnTestSuite3(child)
    
    def lpsnTestSuite4(self, taxO:Taxonomy) -> None:
        """ lpsnTestSuite4:
                Accepts a Taxonomy object as input. Tests if the LPSN data were
                fully amalgamated with the Taxonomy object. Specifically, tests
                that all synonyms were appropriately renamed and that the field
                ncbiNameIsCorrect was accurate. Does not return.
        """
        # constants
        INVALID_STR = "(invalid name)"
        SYN_SUFFIX  = "_syn"

        # no lpsn dictionary for domains
        if taxO.rank < Taxonomy.DOMAIN:
            # get the necessary dictionaries
            rankD:dict = TestTaxonomy.LPSN_D[str(taxO.rank)]
            synD:dict  = TestTaxonomy.LPSN_D[str(taxO.rank) + SYN_SUFFIX]

            # if ncbiName is a synonym, then sciName should have been updated
            if taxO.ncbiName in synD.keys():
                self.assertEqual(taxO.sciName, synD[taxO.ncbiName])
            
            # otherwise, if the ncbiName is correct ...
            elif taxO.ncbiName in rankD.keys():
                # ... then it should match sciName and list the correct boolean
                self.assertTrue(taxO.ncbiNameIsCorrect)
                self.assertEqual(taxO.ncbiName, taxO.sciName)

            # otherwise, it should be marked as invalid
            else:
                self.assertIn(INVALID_STR, taxO.sciName)
        
        # recurse on the children
        for child in taxO.getChildren(set):
            self.lpsnTestSuite4(child)
    
    def pickIngroupTestSuite(self, taxO:Taxonomy) -> None:
        """ pickIngroupTestSuite:
                Accepts a Taxonomy object as input. Tests for the expected fun-
                ctionality of Taxonomy._pickIngroup(). Specifically, confirms
                the function returns a list of species-level Taxonomy objects
                that are internal and have assemblies available for download,
                and that the number of returned objects matches the anticipated
                value.
        """
        # constants
        RANGE = list(range(0, 51))
        MIN_NUM_SEQS = 2

        # copy taxO to check for any modifications
        copyO = taxO.copy()

        # for the range specified
        for numSeqs in RANGE:
            # make sure it fails when fewer than 2 sequences are requested
            if numSeqs < MIN_NUM_SEQS:
                try:
                    # call the function
                    ingroup = taxO._pickIngroup(numSeqs)
                    fail = False
                except: fail = True
            
                self.assertTrue(fail)
            
            else:
                # call the function
                ingroup = taxO._pickIngroup(numSeqs)

                # igroup should be a list
                self.assertIsInstance(ingroup, list)

                # for each item of the list
                ingroupMember:Taxonomy
                seen = set()
                for ingroupMember in ingroup:
                    # must be a Taxonomy object
                    self.assertIsInstance(ingroupMember, Taxonomy)

                    # must be an internal species with an assembly
                    self.assertEqual(ingroupMember.rank, Taxonomy.SPECIES)
                    self.assertFalse(ingroupMember.isExternal)
                    self.assertNotEqual(ingroupMember.assemblyFtp, None)
                    self.assertIsInstance(ingroupMember.assemblyFtp, str)

                    # must not have been seen already (all unique)
                    self.assertNotIn(ingroupMember.taxid, seen)

                    # update the seen set
                    seen.add(ingroupMember.taxid)

                # the list must be the expected size (or empty)
                if ingroup != []:
                    if len(ingroup) < numSeqs:
                        # fewer than requested should return all possible cands
                        self.assertEqual(len(ingroup), \
                                     len(taxO.getAllIngroupCandidateSpecies()))

                    else: self.assertEqual(len(ingroup), numSeqs)
            
            # the calling object should not have been modified
            self.assertEqual(taxO, copyO)

    def pickOutgroupTestSuite(self, taxO:Taxonomy) -> None:
        """ pickOutgroupTestSuite:
                Accepts a Taxonomy object as input. Tests for the expected fun-
                ctionality of Taxonomy._pickOutgroup(). Specifically, confirms
                the function returns a species-level Taxonomy object that also
                has an assembly available for download.
        """
        # constants
        ASSEMBLY_EXT = '.gbff.gz'
        LEN_ASS_EXT  = len(ASSEMBLY_EXT)
        ASSEMBLY_PRE = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/'
        LEN_ASS_PRE  = len(ASSEMBLY_PRE)

        # it's possible pickOutgroup will modify taxO to find an outgroup
        # make a copy so that we can know if a modification occurred
        copyO = taxO.copy()

        # pick an outgroup
        outgroup:Taxonomy
        taxO, outgroup = taxO._pickOutgroup(TestTaxonomy.LPSN_D)

        # the object should still be consistent
        self.assertTrue(taxO._isConsistent())

        # the object should not have '0' as a taxid
        self.assertNotIn(0, taxO)

        # outgroup should never be None
        self.assertNotEqual(outgroup, None)

        # outgroup should be a Taxonomy object at the species level
        self.assertIsInstance(outgroup, Taxonomy)
        self.assertEqual(outgroup.rank, Taxonomy.SPECIES)

        # outgroup has a string for its assembly
        self.assertNotEqual(outgroup.assemblyFtp, None)
        self.assertIsInstance(outgroup.assemblyFtp, str)

        # get the ingroup candidates from before and after calling pickOutgroup
        oriIngroup = copyO.getAllIngroupCandidateSpecies(set)
        newIngroup = taxO.getAllIngroupCandidateSpecies(set)
        
        # the ingroup candidates should remain unchanged
        self.assertEqual(oriIngroup, newIngroup)
      
        # the string should have the expected file extension and prefix
        self.assertEqual(outgroup.assemblyFtp[-LEN_ASS_EXT:], ASSEMBLY_EXT)
        self.assertEqual(outgroup.assemblyFtp[:LEN_ASS_PRE], ASSEMBLY_PRE)
    
        # if taxO has been modified, it should still pass all other tests
        if taxO != copyO:
            self.suiteCaller1(taxO.getRoot())

            # do not recurse on the root! Causes problems with pickOutgroup?
            self.suiteCaller2(taxO)
        
        # pickOutgroup should return the same value if ran a second time
        # taxO should no longer require modification to find an outgroup
        copyO = taxO.copy()
        taxO, newOutgroup = taxO._pickOutgroup(TestTaxonomy.LPSN_D)
        self.assertEqual(outgroup, newOutgroup)
        self.assertEqual(taxO, copyO)

            

    # helper functions for calling the tests in the right order
    def suiteCaller1(self, taxO:Taxonomy) -> None:
        """ suiteCaller1:
                Accepts a Taxonomy object as input. Calls non-recursive test
                suites on the object. Does not return.
        """
        # call non-recursive test suites
        self.equalityCopyTestSuite1(taxO)
        self.equalityCopyTestSuite2(taxO)
        self.equalityCopyTestSuite3(taxO)
        self.equalityCopyTestSuite4(taxO)
        self.loadSaveTestSuite(taxO)
        self.assemblyTestSuite(taxO)
        
    def suiteCaller2(self, taxO:Taxonomy) -> None:
        """ suiteCaller2:
                Accepts a Taxonomy object as input. Calls recursive test suites
                and other tests only suitable for the root of an object. Does
                not return.
        """
        # call recursive test suites
        self.typeMaterialTestSuite(taxO)
        self.lpsnTestSuite1(taxO)
        self.lpsnTestSuite2(taxO)
        self.lpsnTestSuite3(taxO)
        self.lpsnTestSuite4(taxO)

        # if pickOutgroup made taxid a synonym, then update taxO
        # (pickOutgroupTestSuite can call suiteRunner2)
        if taxO.taxid in taxO.taxidSynonyms.keys():
            parent = taxO.parent
            taxO = parent.getDescendantByTaxId(taxO.taxid, resolveSynonym=True)
        
        # call pickIngroup test suite
        self.pickIngroupTestSuite(taxO)

        # call pickOutgroup test suite
        self.pickOutgroupTestSuite(taxO)

    def runner(self, taxO:Taxonomy) -> None:
        """ runner:
                Accepts a Taxonomy object as input. Calls the suiteCaller func-
                tions on the object, which will call an equivalent set of tests
                for all inputs. Does not return.
        """
        # make sure the object was properly constructed
        self.taxonomyConstructionTestSuite(taxO)

        # run the non-recursive tests on the calling object
        self.suiteCaller1(taxO)

        # if taxO is the root, then test one of its children
        root = taxO.getRoot()
        if taxO == root:
            self.suiteCaller1(taxO.getChildren(list)[0])
        
        # if taxO is not the root, then test the root
        else:
            self.suiteCaller1(root)

        # run the recursive tests on the object; do not run on root!
        self.suiteCaller2(taxO)




    # running all the test suites on each test case
    def testA_preBuilt(self) -> None:
        """ Tests a known valid object. Should be identical to testCase5. The 
            only difference is that testCase5 is built from scratch instead of
            read from a table format.
        """
        logging.info(TestTaxonomy.testA_preBuilt.__name__)
        taxO = Taxonomy._loadFromTable(TestTaxonomy.TAX_TABLE)
        self.runner(taxO)
    
    def testB_testCase1(self) -> None:
        """ Tests a Taxonomy object built from Staphylococcus aureus DSM 20231T
            (GCF_001027105.1) blastn results.
        """
        logging.info(TestTaxonomy.testB_testCase1.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_1)
        self.runner(taxO)

    def testC_testCase2(self) -> None:
        """ Tests a Taxonomy object built from Bacillus anthracis str. 'Ames 
            ancestor' (GCF_000008445.1) blastn results.
        """
        logging.info(TestTaxonomy.testC_testCase2.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_2)
        self.runner(taxO)

    def testD_testCase3(self) -> None:
        """ Tests a Taxonomy object built from Ruegeria pomeroyi DSS-3 
            (GCF_000011965.2) blastn results.
        """
        logging.info(TestTaxonomy.testD_testCase3.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_3)
        self.runner(taxO)
    
    def testE_testCase4(self) -> None:
        """ Tests a Taxonomy object built from Methanosarcina mazei S-6
            (GCF_000970205.1) blastn results.
        """
        logging.info(TestTaxonomy.testE_testCase4.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_4)
        self.runner(taxO)

    def testF_testCase5(self) -> None:
        """ Tests a Taxonomy object built from Sulfolobus islandicus Y.N.15.51
            (GCF_000022485.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testF_testCase5.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_5)
        self.runner(taxO)

    def testG_testCase6(self) -> None:
        """ Tests Taxonomy object built from Campylobacter concisus ATCC 33237
            (GCF_001298465.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testG_testCase6.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_6)
        self.runner(taxO)

    def testH_testCase7(self) -> None:
        logging.info(TestTaxonomy.testH_testCase7.__name__)
        """ Tests Taxonomy object built from Salmonella enterica LT2
            (GCF_000006945.2) 16S blastn results.
        """
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_7)
        self.runner(taxO)

    def testI_testCase8(self) -> None:
        """ Tests Taxonomy object built from Geobacter lovleyi SZ
            (GCF_000020385.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testI_testCase8.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_8)
        self.runner(taxO)

    def testJ_testCase9(self) -> None:
        """ Tests Taxonomy object built from Eubacterium limosum ATCC 8486T
            (GCF_000807675.2) 16S blastn results.
        """
        logging.info(TestTaxonomy.testJ_testCase9.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_9)
        self.runner(taxO)

    def testK_testCase10(self) -> None:
        """ Tests Taxonomy object built from Pyrococcus furiosus DSM 3638
            (GCF_008245085.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testK_testCase10.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_10)
        self.runner(taxO)

    def testL_testCase11(self) -> None:
        """ Tests Taxonomy object built from Cyanobacterium stanieri PCC 7202
            (GCA_000317655.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testL_testCase11.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_11)
        self.runner(taxO)

    def testM_testCase12(self) -> None:
        """ Tests Taxonomy object built from Mycobacterium tuberculosis H37Rv
            (GCF_000195955.2) 16S blastn results.
        """
        logging.info(TestTaxonomy.testM_testCase12.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_12)
        self.runner(taxO)
    
    def testN_testCase13(self) -> None:
        """ Tests Taxonomy object built from Methanocaldococcus jannaschii DSM
            2661 (GCF_000091665.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testN_testCase13.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_13)
        self.runner(taxO)
    
    def testO_testCase14(self) -> None:
        """ Tests Taxonomy object built from Mycoplasmoides pneumoniae FH
            (GCF_001272835.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testO_testCase14.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_14)
        self.runner(taxO)
    
    def testP_testCase15(self) -> None:
        """ Tests Taxonomy object built from Sphingomonas palmae JS21-1
            (GCF_900109565.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testP_testCase15.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_15)
        self.runner(taxO)
    
    def testQ_testCase16(self) -> None:
        """ Tests Taxonomy object built from Komagataeibacter hansenii ATCC
            23769 (GCF_000164395.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testQ_testCase16.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_16)
        self.runner(taxO)
    
    def testR_testCase17(self) -> None:
        """ Tests Taxonomy object built from Borreliella burgdorferi B31_NRZ
            (GCF_002151505.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testR_testCase17.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_17)
        self.runner(taxO)
    
    def testS_testCase18(self) -> None:
        """ Tests Taxonomy object built from Williamwhitmania taraxaci A7P-90m
            (GCF_900096565.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testS_testCase18.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_18)
        self.runner(taxO)
    
    def testT_testCase19(self) -> None:
        """ Tests Taxonomy object built from Chlamydia trachomatis A/HAR-13 
            (GCF_000012125.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testT_testCase19.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_19)
        self.runner(taxO)
    
    def testU_testCase20(self) -> None:
        """ Tests Taxonomy object built from Deinococcus radiodurans DSM 20539
            (GCF_000685985.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testU_testCase20.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_20)
        self.runner(taxO)
    
    def testV_testCase21(self) -> None:
        """ Tests Taxonomy object built from Acetomicrobium hydrogeniformans 
            OS1 (GCF_000160455.2) 16S blastn results.
        """
        logging.info(TestTaxonomy.testV_testCase21.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_21)
        self.runner(taxO)

    def testW_testCase22(self) -> None:
        """ Tests Taxonomy object built from Propionibacterium acidifaciens DSM
            21887 (GCF_000426605.1) 16S blastn results.
        """
        logging.info(TestTaxonomy.testW_testCase22.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_22)
        self.runner(taxO)
    
    def testX_testCase23(self) -> None:
        """ Tests a troublesome set of taxids from a phylogenetic marker blast.
            Bug is resolved, but these ids are still a useful test case.
        """
        logging.info(TestTaxonomy.testX_testCase23.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_23)
        self.runner(taxO)

    def testY_testCase24(self) -> None:
        """ Tests a troublesome set of taxids from a phylogenetic marker blast.
            Inconsistent taxonomy after initial construction
        """
        logging.info(TestTaxonomy.testY_testCase24.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_24)
        self.runner(taxO)
    
    def testZ_testCase25(self) -> None:
        """ Tests a troublesome set of taxids identified by a user that caused
            Burkholderia reconciliation to fail during sibling processing.
        """
        logging.info(TestTaxonomy.testZ_testCase25.__name__)
        taxO = TestTaxonomy.buildTaxO(TestTaxonomy.BLAST_TAX_25)
        self.runner(taxO)