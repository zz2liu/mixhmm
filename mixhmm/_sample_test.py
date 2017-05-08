#got from sidcon simu 9 states, each 10 snp
#O, L, LL, LLL, LLLL, LR, LLR, LLLR, LLRR
sample_test = """\
Name	Chr	Loc	BAF	LRR
rs3094315	1	742429	0.00127128240445	-0.544878963982
rs12562034	1	758311	0.510492781795	-0.70610495505
rs3934834	1	995669	1.0	-0.733310545572
rs9442372	1	1008567	1.0	-0.439285042302
rs3737728	1	1011278	0.017996767536	-0.587884846624
rs11260588	1	1011521	0.482992907302	-0.585011900737
rs6687776	1	1020428	0.512747199624	-0.644145737428
rs9651273	1	1021403	0.502184272301	-0.581139386496
rs4970405	1	1038818	0.480399126395	-0.7062926659
rs12726255	1	1039813	0.00170052574618	-0.685544902159
rs11807848	1	1051029	0.683464555045	-0.0651624516101
rs9442373	1	1052501	0.319071386112	-0.252797910048
rs2298217	1	1054842	1.0	-0.394544986169
rs12145826	1	1055892	0.349466014986	-0.197034760367
rs4970357	1	1066927	0.360222251535	-0.312051048571
rs9442380	1	1077546	1.0	-0.175603546455
rs7553429	1	1080420	1.0	-0.270170982137
rs4970362	1	1084601	0.0	-0.332365577278
rs9660710	1	1089205	1.0	-0.247110870798
rs4970420	1	1096336	0.683853521742	-0.314534176309
rs1320565	1	1109721	0.269673482605	0.074159351109
rs11260549	1	1111657	0.23848606275	-0.183890890597
rs9729550	1	1125105	0.746628830732	0.0730821020433
rs11721	1	1142494	0.273945318309	0.209376812023
rs2887286	1	1145994	0.0	0.0922281885836
rs3813199	1	1148140	0.745428663988	-0.106454056604
rs3766186	1	1152298	0.743386773275	-0.0383927884754
rs7515488	1	1153667	1.0	-0.0620714983995
rs715643	1	1162770	0.223216865985	-0.0545604674659
rs6675798	1	1166460	0.761927572242	0.254371611358
rs7524470	1	1182378	0.988007686745	0.309979256441
rs11804831	1	1184667	0.804702599317	0.223858477322
rs6685064	1	1201155	0.0317760988814	0.209247961236
rs3737717	1	1231947	0.180420316546	0.346195639581
rs2765033	1	1300787	0.00655948792223	0.211429533893
rs2649588	1	1303878	0.215911794475	0.158622382418
rs819980	1	1415563	0.203965518107	0.11083679382
rs9439462	1	1452629	0.0	0.062365163022
rs3766178	1	1468043	0.978176603551	0.228496369389
rs2031709	1	1475847	0.00430300807168	-0.0156977203365
rs3128342	1	1476697	0.0147747050893	0.309877473476
rs880051	1	1483590	0.169591779981	0.201075401534
rs2296716	1	1487687	1.0	0.407587505143
rs6603793	1	1495118	0.0	0.22891178652
rs7520996	1	1498897	0.0	0.363195276213
rs3737622	1	1610720	0.83749615842	0.238033372898
rs6603811	1	1695996	0.00643845780379	0.313756343259
rs7531583	1	1696020	0.844743043398	0.3760355685
rs16825336	1	1735586	0.0103138385454	0.394446487852
rs2180311	1	1738594	0.15482265871	0.599025741464
rs6681938	1	1771080	0.502645904349	0.00894142463002
rs10907192	1	1782971	0.999153738313	-0.0164308020107
rs4648592	1	1790894	1.0	-0.0520106769978
rs7525092	1	1799950	0.469720974089	0.0760858935039
rs2474460	1	1833906	0.0	0.235086167552
rs12758705	1	1863485	0.524513910358	-0.211088554505
rs2803329	1	1864441	0.0104741083486	0.299559402412
rs3820011	1	1878053	0.00988151729193	-0.0462685401064
rs2803291	1	1882185	0.0	-0.0509091426212
rs2459994	1	2013924	1.0	0.0157741356943
rs12755035	1	2016221	0.981996850095	0.192396085843
rs884080	1	2016609	0.00185335971313	0.156596953291
rs908742	1	2023116	0.997988130084	0.336565276823
rs4648808	1	2030796	0.00744883921985	0.281653409474
rs6603813	1	2032940	0.400693130907	0.257250260876
rs3107151	1	2041373	0.386424956576	0.282511706082
rs3128291	1	2047883	0.0298002315636	0.20306840321
rs3128296	1	2058766	0.375468673338	0.140296822492
rs3753242	1	2059541	0.424579264835	0.205138978061
rs424079	1	2061200	0.594882294987	0.232360760431
rs262669	1	2072349	1.0	0.154316310089
rs2257182	1	2072426	0.0	0.222299744358
rs262654	1	2079386	1.0	0.278629867918
rs3052	1	2086498	0.333968046264	0.166885734519
rs262688	1	2103425	0.665069635065	0.418303122272
rs2460002	1	2109693	0.339719087834	0.24335190803
rs6665593	1	2130121	1.0	0.249281371288
rs7512482	1	2136826	0.0184950947708	0.353986293422
rs262683	1	2145729	1.0	0.309569983883
rs2460000	1	2146222	0.0174727684167	0.380347356653
rs263526	1	2163364	0.515204342619	0.49214393948
rs10910034	1	2165898	0.0	0.208111460235
rs1713712	1	2166021	0.500866418636	0.21060642128
rs260513	1	2170384	0.491761686239	0.431450330691
rs260512	1	2172330	0.0	0.403374714434
rs10797417	1	2182293	0.0	0.405416951023
rs7547453	1	2184977	0.507190484393	0.417387939644
rs11588312	1	2188917	0.486765987838	0.404182326278
rs6673129	1	2192634	0.503095222896	0.449916432423
rs7553178	1	2194615	0.503929829014	0.363519479579
"""

