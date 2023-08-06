classdef sos < handle
   % sum-of-sinuoids autocorrelation model
   properties
       name
       corr_dist          % the correlation distance
       sos_wn             % the angle frequencies of sinusoids
       sos_an             % the amplitude of sinusoids
       sos_phase          % the random phase [-pi, pi]
   end
   properties(Constant)
       def_corr_dist = 10 % the default correlation distance
   end
   properties(Dependent)
       Num_sin            % the number of the sinusoids
       Dim                % dimessions
   end
   
   methods
       function obj = sos(varargin)
           if isempty(varargin)
               obj.getDefault();
           else
               obj.getParameters(varargin{1});
               if length(varargin) == 2
                   obj.sos_phase = 2*pi*(rand(obj.Num_sin,varargin{2})-0.5);
               end
           end
       end
       
       function getDefault(obj)
           obj.sos_an = [0.0034051405,0.0041435524,0.0031057296,0.0037385155,0.0028595286,0.0035250664,0.0042655719,0.0041326880,0.0033275280,0.0026985114,0.0031031505,0.0032311175,0.0031204903,0.0035986577,0.0033520667,0.0031245584,0.0032015278,0.0028551649,0.0029928710,0.0032212916,0.0031765285,0.0030485447,0.0029924137,0.0033869776,0.0035793406,0.0033068529,0.0029074550,0.0031569414,0.0038354655,0.0034465210,0.0037632817,0.0034626327,0.0031276117,0.0030877416,0.0044221752,0.0033471494,0.0028564685,0.0039819521,0.0036830767,0.0030967339,0.0037964436,0.0030617260,0.0035739061,0.0028207228,0.0027614194,0.0027898883,0.0032352151,0.0032510718,0.0028670961,0.0035539882,0.0028287955,0.0030324569,0.0035560899,0.0032248734,0.0028305505,0.0035705706,0.0028080032,0.0037677353,0.0029637106,0.0034508877,0.0033473717,0.0038659100,0.0031471662,0.0034851443,0.0029131453,0.0033801424,0.0035844899,0.0035273228,0.0036552723,0.0035665664,0.0038018306,0.0027248212,0.0037908757,0.0033208595,0.0028762184,0.0030066122,0.0029976671,0.0028561582,0.0035096908,0.0031469685,0.0030656010,0.0028405162,0.0036586691,0.0027223479,0.0035046346,0.0029932936,0.0036951040,0.0030293223,0.0029946133,0.0038888948,0.0036677381,0.0035332816,0.0031438575,0.0035028346,0.0036068810,0.0030046566,0.0035694079,0.0035652376,0.0032016782,0.0034914839,0.0033734213,0.0035703678,0.0037424434,0.0027377631,0.0030122222,0.0040011997,0.0030659870,0.0029056582,0.0036222022,0.0034474749,0.0030788868,0.0035164992,0.0030350373,0.0035398128,0.0032384330,0.0029845366,0.0034313311,0.0032696007,0.0025769221,0.0028864718,0.0026076518,0.0030266980,0.0039242203,0.0031996693,0.0031233779,0.0032771337,0.0036989853,0.0029617595,0.0034855823,0.0034916764,0.0028771944,0.0034192097,0.0031328814,0.0028454650,0.0031736335,0.0037752874,0.0036541242,0.0034098320,0.0031306997,0.0032155802,0.0035315705,0.0034346541,0.0031123280,0.0028599098,0.0036730391,0.0034071642,0.0029522751,0.0035981741,0.0028447176,0.0035065264,0.0030230850,0.0034797883,0.0032159064,0.0035631417,0.0028259994,0.0032114936,0.0037530577,0.0034895046,0.0035094509,0.0029437637,0.0033659504,0.0030139375,0.0031993436,0.0029800392,0.0032246520,0.0039470582,0.0039869365,0.0033815857,0.0038972902,0.0029015162,0.0039943429,0.0035512266,0.0035539484,0.0034451897,0.0035742037,0.0031548666,0.0027990350,0.0032916607,0.0035920097,0.0027901221,0.0036838267,0.0029560227,0.0035254972,0.0037635209,0.0030795815,0.0032856374,0.0033707737,0.0028907803,0.0030574079,0.0031704786,0.0029156916,0.0029763754,0.0029605948,0.0033430008,0.0032291978,0.0027545465,0.0033353711,0.0031965026,0.0035449099,0.0037247760,0.0032931725,0.0038662776,0.0029822432,0.0033817866,0.0029727868,0.0035547770,0.0037200793,0.0029865792,0.0033926831,0.0032877557,0.0037862109,0.0029375437,0.0027018311,0.0033585096,0.0034984224,0.0033204819,0.0034381056,0.0034457115,0.0032334919,0.0035735797,0.0039901906,0.0030845772,0.0036994203,0.0034074713,0.0030894203,0.0033219485,0.0034124909,0.0034935861,0.0029881559,0.0035537842,0.0036224083,0.0035811558,0.0041373563,0.0029412173,0.0035820066,0.0033957732,0.0035834610,0.0036618568,0.0037487266,0.0039420687,0.0033655171,0.0038635924,0.0029204516,0.0031267796,0.0036796900,0.0038633896,0.0037311087,0.0035774999,0.0034190898,0.0030454835,0.0032051082,0.0035514652,0.0035769427,0.0032843875,0.0033552218,0.0035816976,0.0038746542,0.0030809727,0.0032694694,0.0031441653,0.0031272708,0.0031805558,0.0029248982,0.0034937290,0.0031226755,0.0028838629,0.0030040492,0.0037344822,0.0035046099,0.0032657592,0.0032463525,0.0031584443,0.0033203682,0.0036609240,0.0037756029,0.0036562260,0.0035896003,0.0034311439,0.0033922698,0.0034424225,0.0029122678,0.0037012435,0.0030703486,0.0028418843,0.0035933936,0.0035824375,0.0035293996,0.0034632178,0.0032013380,0.0036675781,0.0032025177,0.0033752944,0.0030456551,0.0030675537,0.0035237100,0.0030882352,0.0035816713,0.0035344122,0.0034360606,0.0037708699];
           obj.sos_an = sqrt(2*obj.sos_an);
           obj.sos_wn = [-0.064455755,0.013439692,-0.056105554;-0.21950509,0.096034914,-0.25416169;-0.031303219,0.13855690,-0.028452408;-0.14131607,-0.019121891,-0.22707486;0.019583553,-0.10127676,-0.083872974;0.061835412,0.016713411,-0.020435937;0.096401416,-0.39010590,0.18863614;0.012868898,0.046033647,-0.31420416;0.24783459,0.17127301,0.13191184;-0.081635170,-0.089179941,0.091116555;0.047289278,-0.10471740,0.024800647;0.24610317,0.076707400,0.11226680;0.034956854,0.031548928,-0.095792159;-0.019295162,-0.013324788,0.030870482;0.26522374,0.12083730,-0.19132788;0.055366885,0.091838226,-0.024166400;0.12484374,0.12240096,-0.15754502;0.043094188,0.087376006,-0.080753379;0.068010069,0.016318118,-0.089239225;0.11558214,-0.16440074,0.091007397;0.063589349,-0.12966773,0.094766669;0.038268596,-0.11140698,-0.028635757;0.065409541,-0.10985204,-0.14690399;0.21951044,0.075641476,-0.14845857;-0.24037379,-0.28116682,-0.26402804;-0.18282151,0.051575124,0.036800697;-0.14471464,-0.0063504176,0.037104499;0.0070047211,-0.081930228,0.064997494;-0.022847215,-0.32178074,-0.15817790;0.085803621,0.16897658,0.34743434;-0.23439683,0.16408034,-0.17936440;-0.24672596,0.082264967,0.12003037;-0.21263404,-0.041185938,-0.054952700;0.048091698,-0.088524014,-0.051276818;0.12414955,-0.44222212,-0.081647053;0.070179276,-0.052450854,0.041166432;-0.036763150,0.023417221,-0.14887901;0.044956900,-0.21667616,0.22807589;-0.26999411,-0.072571389,0.047277860;0.053148065,-0.046734296,-0.087295435;-0.35882023,0.010249967,0.15064116;-0.076142989,0.069360666,0.055397805;0.0099386470,-0.012644349,0.056663696;0.12481805,0.13362192,-0.081831746;0.13100091,-0.036292061,-0.10364211;-0.031384677,-0.00048892590,-0.13084799;-0.32823354,0.16678141,-0.25736496;-0.16803236,-0.27692837,0.21741344;0.10946541,0.013008613,0.12162443;-0.0036848858,-0.039140530,-0.046528142;0.036750458,0.050658345,0.12882870;0.059557360,-0.066955686,-0.075388253;-0.17983298,0.23123744,-0.063930005;-0.058951423,-0.19299892,-0.00046347050;-0.027432406,0.050130580,0.45108721;-0.037368666,0.018971298,0.048908420;-0.0095554953,0.036676954,-0.14251523;-0.21434906,0.21722488,0.20001404;0.046888344,0.089361385,0.072574370;-0.11944297,-0.22210455,-0.041989923;0.081054777,-0.035079576,0.041341845;-0.31631291,0.22552937,-0.0086688325;0.16184682,0.14158405,0.036902294;0.19263659,0.25597450,-0.0048904433;-0.046801217,-0.081074044,0.15816684;-0.0083855530,0.081031442,-0.038897436;-0.013864271,0.024874195,0.049357105;0.069010012,-0.011994654,-0.0061423155;-0.00022581631,0.33381596,0.065264352;-0.012092848,-0.033365373,0.0018169039;-0.042661514,-0.41800219,0.011240603;0.10337584,0.053941902,-0.086158484;-0.38693455,-0.20289235,-0.00010028003;-0.15679963,0.15956497,-0.046156142;0.064687170,0.15128826,-0.039175101;0.058925938,0.18220654,-0.044759430;-0.10364697,0.090353720,-0.067607656;0.043268915,0.13685974,0.099037066;0.35781318,-0.033824753,0.15794691;-0.16235340,0.083577596,0.095003501;-0.10745410,0.0026796903,0.19167905;0.0080160815,0.064146228,-0.12755299;0.025037047,-0.099349245,-0.27747223;0.0015462569,-0.0064206738,0.13922441;-0.069976777,-0.058792721,-0.21751878;0.15129593,0.066551276,0.030167188;0.074662797,0.37468910,0.26904181;-0.015338658,0.10543114,0.051341254;-0.062027790,-0.10702091,-0.032238714;-0.39807832,0.041407671,0.077866025;-0.11478776,0.31000063,0.33808172;-0.049236283,0.028801057,0.045206785;0.0017183678,0.048507467,-0.093988508;-0.0050686132,0.035604481,0.066079699;0.13705559,-0.078986414,0.14949496;0.16007648,0.050429970,0.032201461;-0.11392801,0.20466323,0.061685786;0.050633471,0.0019579113,-0.0040332549;-0.058041874,-0.14330927,-0.41737655;-0.29586947,-0.0078383470,-0.014332021;0.068456411,-0.17820145,-0.0063378620;0.019859599,0.023117239,0.040387787;0.18437976,-0.23832361,-0.075346351;0.059387393,0.12924396,-0.069162950;-0.039033223,0.098876469,0.053491902;0.19378607,-0.010251210,0.17855862;-0.035457127,-0.18259120,-0.088837959;0.084596708,-0.045417741,-0.17742181;0.14010000,0.14357036,-0.34585384;0.047884136,-0.19075248,-0.0025211961;0.022178944,0.21996053,-0.096612126;0.12491512,-0.038810771,0.23607524;-0.036472242,-0.18509951,0.042478807;-0.33381683,-0.082636118,-0.13892443;-0.021019913,-0.10111051,-0.026189514;0.039849650,-0.14178605,-0.072870612;-0.0053715417,0.040164668,0.075057775;-0.029063996,-0.097179942,-0.0033870458;0.059124142,-0.032899000,-0.12645970;0.099136166,-0.074658461,-0.068439454;-0.098569326,0.020772882,0.11819197;0.13169020,-0.00046980241,0.056255907;-0.13232449,-0.25809839,0.30636171;0.098758288,0.19013956,0.14350122;-0.057497341,0.15405791,-0.070611492;-0.19654471,0.036511067,0.099517800;-0.28634039,-0.19826363,-0.069728836;0.0043304050,-0.072709076,0.088505372;0.067156754,0.0062551619,0.031831548;-0.015916759,0.0021009336,-0.074294873;0.022427974,-0.054751672,0.12212449;-0.16263694,0.15458684,0.057057358;-0.091703244,0.099124417,0.0051252372;0.039032813,-0.10976655,-0.096256442;-0.058679581,0.026078247,-0.085187912;0.0083274031,0.21221358,0.20943788;0.040845565,-0.30014476,0.022321703;-0.18563096,0.079414897,0.040159736;0.12531917,0.0015322963,0.021106206;-0.046923857,0.066363715,-0.063777320;0.065955475,-0.018629240,-0.021382989;0.27650189,-0.0028397394,-0.20245983;0.21203077,-0.018805647,0.059052728;-0.031761318,-0.10327688,-0.11998299;0.029812859,0.34285250,-0.074710272;-0.20032614,0.20253961,-0.30787769;-0.14178689,-0.047112606,0.0076201870;0.032200161,0.0021708913,-0.028293969;-0.041822996,0.041318402,-0.11486190;0.048625082,-0.047129869,0.030240489;0.10739152,-0.074411280,-0.039132282;-0.038887754,-0.056446381,0.018892178;0.16375220,-0.099353030,0.048352066;-0.010257031,-0.038837016,-0.018267632;-0.16463006,0.16934520,0.41575018;0.092129201,-0.18610433,-0.14101399;0.17972815,0.27550519,0.14403038;0.33356011,0.040010624,-0.14453723;-0.23831160,0.28121409,-0.15254299;0.090713136,-0.014729225,-0.075562932;-0.0016411705,0.0072443453,0.091301486;0.10877039,-0.059376672,-0.037375499;0.22989672,0.055963945,0.0051861140;0.12804514,0.056122590,-0.0067978818;-0.031816259,0.22999687,-0.097380631;0.27841291,0.0091954758,0.29371056;0.22443864,-0.37060821,0.080176741;0.062853180,0.23305202,-0.15076706;0.058413472,0.23009326,-0.22990108;-0.051455341,-0.10558267,-0.12511301;-0.15613225,-0.20132077,0.25935927;-0.0064839046,0.051445868,-0.0012007713;0.00026647258,0.047394574,-0.037232652;-0.11524726,2.0747360e-05,0.30389947;-0.011884146,0.031037146,-0.0073397039;0.12057044,0.014353521,0.019784456;-0.0093778186,0.036009744,-0.14079139;0.064473651,-0.053935021,0.16892064;-0.010016029,-0.013405675,-0.037003927;-0.12084470,-0.086178824,0.068026140;0.25049961,-0.21464613,0.012864460;0.072224975,-0.0041677658,-0.089953005;0.070732355,-0.11048662,0.13776553;-0.28933620,-0.30398095,0.076618522;0.14265122,0.13486932,0.13146265;0.12076825,0.20716372,-0.057086285;-0.073030733,-0.029566932,0.044475950;-0.067417584,-0.029524546,-0.14636970;0.13142310,-0.037158240,-0.0051325285;0.14918767,0.15810665,0.037214290;-0.010082254,0.077573515,0.097737774;0.058368150,0.15051414,0.040586714;0.18102108,-0.033503316,0.081750646;0.013177372,-0.042172704,0.081980012;-0.38623720,-0.00043652969,-0.25394568;-0.045988839,-0.13867836,0.069959335;0.083916202,-0.052236084,-0.0071570431;-0.19527639,0.067331679,-0.055327091;0.032223541,0.26162758,0.048866875;0.12464774,-0.17093414,0.20802990;-0.047411982,-0.085621156,-0.0025438557;-0.31753299,0.094216660,-0.084396116;-0.056107771,-0.048590191,-0.091767080;-0.18370169,0.14495331,-0.052754205;0.10683133,-0.013181080,-0.057502996;0.10474759,-0.15123974,-0.17973939;0.050584290,-0.13813940,0.15266404;-0.097285397,0.031208931,-0.079914875;0.059569556,-0.054627907,-0.20577610;-0.028385585,-0.084398843,0.22280473;0.39389238,-0.073751330,0.094793245;0.083746746,-0.11603183,-0.074883096;-0.11529025,-0.081850067,0.096025094;0.050242666,0.074917600,0.022466777;0.16343451,-0.11128823,-0.20238830;-0.065529995,-0.046208028,0.19490327;0.063544460,0.27581668,-0.12525460;0.11566138,2.8426137e-05,-0.30492389;-0.041653592,0.14172661,0.45199209;0.0011364452,0.029854046,0.0054810122;-0.27365583,0.10582364,0.21038282;-0.090637147,-0.048646133,0.053779796;0.24153139,-0.15208137,-0.12566669;-0.16379841,0.15243341,-0.13886385;0.089142784,0.071816742,-0.028138526;-0.19978252,0.090623409,-0.056955278;-0.057991952,-0.13739565,0.24097636;0.11209220,-0.22232030,0.11146279;0.14726582,0.084310815,0.030381316;0.041562721,-0.014541662,-0.052545853;-0.042809401,-0.0040492988,0.38754618;-0.12596358,0.12913261,0.31082541;-0.059559211,0.15788573,0.26936564;-0.14606225,-0.014788949,-0.095853984;-0.0072703147,0.0081958342,0.0024456107;0.22073804,0.10074756,-0.10606863;-0.0023663111,-0.0070356224,0.0089164032;-0.038049549,0.30506444,0.069760613;-0.25514671,0.28435472,0.17413197;-0.19815190,-0.10634181,-0.22590522;0.059360154,0.066942707,1.1979302e-05;0.10721684,0.35535669,0.16854863;0.0076197898,-0.13854238,-0.10555340;0.035697334,0.060012694,0.084455848;-0.27057102,-0.036714412,0.045484446;-0.21789187,-0.11027600,0.24629572;0.33362019,0.18012559,0.25659287;-0.095386565,-0.46583995,0.0027050718;0.089468949,-0.065175772,-0.30865541;-0.087475821,-0.040935289,-0.072100900;0.19600008,-0.32125956,0.27572620;0.038913708,0.028886074,-0.0028635918;-0.010501595,0.034366231,-0.026264047;-0.011121782,0.26797301,-0.063050866;0.034125306,-0.15195790,-0.17035826;-0.00082769943,-0.0033051453,0.0027618890;-0.18291689,-0.012672906,-0.32870260;-0.067258187,-0.0057250764,0.085545875;-0.041346818,-0.10156199,-0.28666276;-0.16978306,-0.14863074,0.011484682;0.042362679,0.084825560,-0.049866132;-0.19814284,-0.044414878,-0.12417739;-0.10010629,-0.081172049,-0.085877366;-0.35220227,-0.095539786,0.15708779;-0.071565546,0.097323552,0.0036949387;0.086769074,0.11085805,0.060751576;0.011475728,0.16155991,0.11776412;0.091872148,0.075555012,-0.34479010;-0.22487818,-0.35014102,0.16819452;0.010542105,-0.20774110,-0.10650031;-0.078833655,-0.056214109,-0.044308677;-0.26615050,-0.24771047,0.15345307;0.078516260,-0.045254052,0.042765807;-0.42688531,0.18590792,0.0068131983;-0.13768072,0.25065210,0.19309506;0.13715616,0.26559138,0.12930423;-0.017263783,-0.017548811,0.022164539;-0.063102439,-0.13474807,-0.26286104;0.11172382,0.052961908,-0.20019868;0.042908989,-0.072776780,-0.0011206710;-0.12460533,-0.066202700,0.035724804;-0.21270381,-0.084235966,0.40021065;0.060256332,-0.098571159,-0.029959869;-0.010677376,-0.081705950,0.097684503;0.19222283,-0.39165661,0.011625860;-0.017894810,-0.0073840693,0.0064534084;-0.12270124,0.037206117,-0.23073874;0.24551915,0.20437759,0.13887672;0.11911546,-0.0014879243,-0.20880902;-0.25041476,-0.18172035,0.061626218;-0.026178543,-0.094074875,0.25384754;0.34008145,-0.11745042,0.24712971;0.017479269,-0.12229910,0.025050597;-0.13359259,-0.039401423,-0.14114435;-0.046961702,0.034831356,0.045681909;0.17638089,-0.027808350,0.0057007140;-0.00027539564,-0.00025916906,-0.0018889387;0.023544936,0.011158262,0.060647037;-0.065736681,-0.035296887,-0.035650194;0.050711870,-0.29617795,0.19335654];
           obj.name   = 'comb_300';
       end
       function getParameters(obj,type)
           switch type
               case 'exp_300'
                   obj.sos_an = [0.0023430123,0.0019328272,0.0054924460,0.0030965735,0.0024775655,0.0059098322,0.0044421386,0.0041632969,0.0032970456,0.0029205927,0.0049209846,0.0033964361,0.0018166252,0.0021206713,0.0021015496,0.0021288141,0.0031100498,0.00092768826,0.0016140407,0.0021448710,0.0020131273,0.0025845580,0.0042274361,0.0025769144,0.0032260800,0.0019239506,0.0048219375,0.0019902387,0.0051851431,0.0071554920,0.0021484247,0.0026543166,0.0015075427,0.0018934862,0.0026696555,0.0019859462,0.0020473329,0.0051634563,0.0050963308,0.0047447849,0.0037167836,0.0036697679,0.0019890731,0.0026601909,0.0047302423,0.0029949665,0.0020576871,0.0023197385,0.0018179885,0.0019994276,0.0031025780,0.0020075236,0.0013988281,0.0038664208,0.0056430516,0.0017406073,0.0017139619,0.0039725443,0.0018365794,0.0050236736,0.0049369684,0.0021155721,0.0020625677,0.0050697997,0.0058309995,0.0023728679,0.0022468893,0.0052874056,0.0052872915,0.0045084530,0.0038450048,0.0048523964,0.0015831288,0.0037139477,0.0029664480,0.0050164675,0.0014916159,0.0033445400,0.0052565355,0.0025272192,0.0043060714,0.0049264636,0.0037795752,0.0038526638,0.0029141251,0.0021053895,0.0062903748,0.0050163944,0.0050615058,0.0049189762,0.0019876685,0.0027095636,0.0053811250,0.0049225092,0.0042002029,0.0024296290,0.0017714619,0.0060036969,0.0038651098,0.0018590153,0.0016007727,0.0019335349,0.0025749269,0.0025326877,0.0041282224,0.0018072544,0.0056832978,0.0025229319,0.0021499649,0.0030743612,0.0045368467,0.0028688838,0.0056863781,0.0014056339,0.0030582158,0.0018763002,0.0031735904,0.0022903928,0.0016099615,0.0023370441,0.0049092928,0.0018406535,0.0037414413,0.0051725670,0.0053524193,0.0024852143,0.0063182376,0.0054864003,0.0026407000,0.0037034166,0.0022725773,0.0016158567,0.0024946006,0.0021480606,0.0034985791,0.0046482771,0.0017057254,0.0016425904,0.0050490703,0.0012893331,0.0042672455,0.0039980658,0.0062505710,0.0036732072,0.0049424525,0.0020902848,0.0022932962,0.0024978057,0.0041638589,0.0015450990,0.0043998999,0.0025410959,0.0044458746,0.0022437985,0.0024071566,0.0015681593,0.0018560443,0.0041967947,0.0049400800,0.0052059088,0.0022037602,0.0055585573,0.0060729096,0.0021506050,0.0025783246,0.0017008084,0.0025540362,0.0012702739,0.0019591500,0.0019441649,0.0020777816,0.0026061093,0.0026597746,0.0020057389,0.0018291180,0.0023936960,0.0045439228,0.0015409937,0.0051995507,0.0017742509,0.0044042338,0.0032576057,0.0039865253,0.0047976281,0.0016043525,0.0047411220,0.0020410214,0.0026327078,0.0017377924,0.0038267884,0.0066146664,0.0043138918,0.0017185615,0.0021834471,0.0040194453,0.0048045223,0.0025297843,0.0024143993,0.0031976446,0.0069940537,0.0015288444,0.0040157046,0.0024925193,0.0041644210,0.0052811047,0.0062476313,0.0019902803,0.0020909249,0.0022577066,0.0044005397,0.0029194478,0.0025982161,0.0049952543,0.0044620200,0.0052444539,0.0054373657,0.0048722364,0.0039962498,0.0038684241,0.0025876591,0.0016908495,0.0048146998,0.00075390312,0.0023070844,0.0045266482,0.0028742638,0.0027938127,0.0019718611,0.0031967340,0.0053316951,0.0018358561,0.0034948380,0.0050771777,0.0020258906,0.0022157016,0.0020375613,0.0024015112,0.0049135336,0.0026509396,0.0071082525,0.0028151248,0.0022877965,0.0041859304,0.0050274255,0.0025370168,0.0024314604,0.0054452643,0.0021641124,0.0024339529,0.0024269754,0.0061778361,0.0019735002,0.0046184319,0.0043608085,0.0020048283,0.0063209124,0.0048462152,0.0022035106,0.0039232029,0.0017139048,0.0021082878,0.0022387463,0.0029105558,0.0045019072,0.0053022713,0.0056340029,0.0052915164,0.0021692009,0.0042305370,0.0025574078,0.0045293798,0.0023174023,0.0042560897,0.0018380883,0.0033342482,0.0012098446,0.0017697137,0.0045352075,0.0033470530,0.0016536782,0.0051572388,0.0017499470,0.0028167274,0.0025729600,0.0017094977,0.0026403512,0.0032702442,0.0023159662,0.0042253495,0.0010261267,0.0040574078,0.0060338131,0.0031104868,0.0032074475,0.0015831080,0.0060359202,0.0035305421,0.0036965120,0.0044913348,0.0039727888];
                   obj.sos_an = sqrt(2*obj.sos_an);
                   obj.sos_wn = [0.043350309,-0.022881916,0.038323190;-0.019704791,0.086650424,0.018369351;-0.68894714,0.26792961,0.49441665;0.077314667,-0.062018644,0.20145677;-0.073889382,0.0055997288,0.053045362;-0.31390917,-0.22297275,-0.67040473;-0.47599781,0.25801206,-0.57126540;-0.029896110,-0.024287788,-0.11841456;0.20950301,0.15686476,-0.51263738;-0.045905106,-0.18226901,0.059103273;0.12237323,0.093140453,-0.095575102;-0.036570575,-0.093996130,0.17182663;0.055196933,-0.064354494,0.040724561;-0.065451950,-0.050741665,0.018112289;-0.00080201350,-0.073176108,0.011520850;-0.062658250,0.022632837,-0.081736207;0.023457287,0.15499812,-0.025944278;-0.24885981,0.22686175,-0.16708992;-0.034577101,-0.047952164,-0.050090399;0.026419882,-0.062982060,0.045774281;0.050411873,-0.11751931,0.079181939;0.0075605637,0.013936137,-0.015845956;0.10901345,0.0010091836,-0.043403104;-5.7441350e-05,0.031507507,0.0063721333;-0.045212436,-0.16402631,0.12632360;0.089230515,-0.048513755,0.032594066;-0.19832383,-0.19980176,-0.066176407;-0.019443793,-0.063731074,0.023920096;-0.19490725,0.018820493,0.095567130;0.73821771,0.32199371,0.50952441;-0.016387526,0.069576934,-0.0096984366;0.0050490620,-0.0066993656,-0.00022946707;-0.10515825,0.42500594,-0.29332900;-0.014974530,-0.0029363709,-0.087938368;-0.022021407,-0.0030532184,0.098576859;-0.0019610797,-0.017289767,0.074515283;-0.045606688,-0.061829969,-0.021046774;0.24145380,0.12109397,-0.18871154;-0.14736566,0.063207112,0.037309475;0.024948834,-0.18828699,-0.053053543;0.038786538,0.23523679,-0.034981698;0.044505585,0.079800971,0.049716048;-0.0030771946,0.010867951,-0.093341269;0.0024464030,-0.00095787155,0.0054412535;-0.15803568,-0.25457272,-0.17323270;-0.10336232,0.15643625,-0.038879301;0.020650705,-0.046856876,0.10500886;-0.082826018,-0.014378116,-0.074470788;0.040180400,0.057155956,0.0042917165;-0.11550473,0.17116641,-0.19965997;-0.069974676,-0.13210247,0.26161614;0.011155738,-0.12317985,0.072927333;-0.054262437,0.031322826,0.048016023;0.32817551,0.50305915,0.25146276;-0.25508559,0.58908015,0.45590830;0.041209403,-0.091772668,0.028560122;-0.060999375,0.063883767,-0.017440883;0.21530290,-0.34554273,0.020850955;-0.0053435583,-0.077553272,0.031850155;0.071181379,-0.047073148,-0.089183360;0.39982834,0.062917292,-0.11529570;-0.023202572,0.12473126,-0.10072944;0.039895535,0.040701810,-0.076053433;-0.22260208,0.66077685,0.14688139;-0.081062153,0.16413344,0.38236597;-0.026299965,-0.014969685,0.018408787;0.027349873,-0.015278524,-0.032673933;-0.036935221,-0.070730731,-0.19244787;0.34514999,-0.19525410,0.88629031;-0.31797701,-0.039844889,-0.080630302;-0.094015926,-0.023096820,0.044285335;0.048857562,0.037417665,0.15053250;-0.070512787,0.00024149028,0.0034085806;0.064300336,0.18048941,-0.18373090;0.13784127,-0.13373640,0.052507248;-0.13822061,-0.027376899,-0.021041604;-0.024765586,0.074906148,0.040445685;0.079694428,0.33894342,-0.020666298;-0.071015425,0.11030377,0.11622803;0.29617307,0.35794407,-0.36317056;-0.14890355,-0.095431313,-0.097145177;-0.14805971,0.015120807,0.54956383;0.17745133,-0.16982944,0.013268190;-0.17018937,0.036053747,-0.065933585;-0.13802390,0.074411683,-0.030168107;0.053298164,-0.027500005,0.074103974;-0.033211499,0.85971326,0.45718497;-0.049838804,-0.11096874,-0.11786230;0.11529204,0.13965544,0.077346906;0.14053325,-0.038112227,0.41376007;-0.032207903,0.058688499,-0.070564047;-0.0011913246,0.12054750,-0.13512996;0.29215634,0.15051423,0.20101452;0.0063891225,0.038250227,-0.29806921;-0.13965704,0.023036120,-0.14665250;0.024387987,0.0097603872,-0.016402548;-0.075318001,0.065361589,-0.028001659;-0.48228562,-0.19993879,0.036986023;-0.45059109,-0.55319434,-0.13816383;-0.021807006,-0.066615880,-0.041271053;0.083857536,-0.19734228,0.35506222;-0.018955482,0.066288173,0.020841027;-0.019401496,0.023170175,-0.0086196447;-0.0043453188,0.0042687468,-0.036625758;-0.23264615,0.11696841,0.023296574;0.0076163439,-0.065654330,-0.078870736;0.024814516,-0.21285723,-0.16849262;0.51593447,-0.28008953,0.23965430;0.042897224,-0.11393737,0.048574720;0.14371522,-0.061575569,0.031070407;0.0071442751,0.13397942,0.039582536;0.0058103595,0.38043210,-0.54355139;-0.14761685,-0.096082114,0.034871738;0.0012857651,-0.067470670,-0.058852553;-0.018241610,-0.049179800,0.13613582;-0.033446893,-0.055958509,0.024028070;0.078908831,0.020344416,0.093692504;0.037275799,-0.42965931,0.45821911;-0.057199873,0.064781807,0.040183809;0.086353846,-0.061837710,-0.021993678;0.20010231,0.27459669,0.41262317;-0.022974528,0.032058354,0.054989889;-0.017409852,0.16175126,-0.0037024082;0.19623207,0.19269951,-0.78617674;0.15755621,0.022213621,-0.024854602;-0.027112499,0.14036891,-0.057880189;-0.092119366,0.29507563,0.45311886;-0.39206797,-0.13513131,-0.24386109;-0.0076723895,-0.010618970,-0.0066918354;-0.32438424,-0.21370755,0.64802307;-0.025975373,-0.045853581,-0.0067560826;0.23319499,-0.18615864,0.10139663;0.015062222,-0.035509296,-0.0072794794;-0.0017912035,0.049413402,-0.049568322;-0.061565354,0.024228116,-0.13363995;-0.17200479,-0.085112788,-0.14940676;-0.064627156,-0.011252055,-0.036244728;0.035221297,0.038482487,0.062509976;-0.079094097,-0.82895827,-0.10349658;0.092040986,0.40310329,-0.30292729;-0.027831076,-0.074737668,-0.32498825;0.59406751,-0.21160364,0.45537826;-0.11651430,-0.56876904,-0.63981086;-0.20923738,0.088334434,-0.39413932;0.11198165,0.058357596,0.0041896543;-0.0062932014,0.026696013,0.053926021;0.024326472,0.019318353,0.049674287;0.10628553,0.0014117782,0.0054714945;-0.069590792,0.0017556213,-0.26521388;0.029651152,0.043981686,0.063904129;0.056937542,0.31653851,0.12842715;-0.040904935,0.043376550,-0.11667871;-0.090316609,-0.084361359,-0.037572663;-0.037550781,0.044241805,-0.0042975629;-0.036104269,0.0091947401,-0.0033893958;0.037001200,-0.0095530804,-0.069402382;-0.051904213,-0.023670360,-0.029545534;-0.086774848,-0.032816347,0.19173677;-0.20543154,-0.055619910,0.29447547;-0.76096648,-0.059166491,-0.27197766;0.024413832,0.11349352,-0.094114028;0.25998163,-0.0045417296,-0.25762755;-0.61314136,-0.37976658,0.015492622;0.035875041,-0.094071686,0.14423344;0.041727860,-0.045521945,0.11954781;-0.047733206,0.082682237,-0.048235703;0.024301570,-0.078830384,-0.070947066;-0.035205774,0.061638135,0.053468831;-0.0050549451,-0.065688983,0.10936599;-0.051411591,0.058811221,-0.061777540;-0.029055027,-0.020748748,0.039757133;0.11629423,-0.022621291,0.010972959;-0.00020845205,-0.00080535171,-0.0046173595;0.054323893,-0.018625541,-0.0077526830;-0.062337127,0.25383380,-0.16588065;0.031010719,0.072851613,0.019323861;0.26159000,-0.23576443,0.015491411;0.045550752,-0.075859584,-0.022409806;-0.042265549,-0.12828970,-0.086748309;0.034226622,-0.072826587,-0.062639721;0.016157299,-0.13326879,-0.044826627;0.041471280,0.12151290,-0.017031696;0.33744740,0.45654985,-0.35494912;0.083207577,0.0082182884,-0.18238567;0.068985872,0.034132160,0.023133958;-0.023089884,-0.25747213,-0.033562593;-0.0058802776,0.076461129,-0.012981180;-0.0048665451,0.011813949,0.0078727510;0.074752025,-0.018081976,0.018785609;-0.18776064,-0.53383780,-0.033222511;-0.22138591,0.42391458,0.88536149;-0.42982915,-0.58592665,-0.43341917;0.0078118797,-0.031921864,-0.069042757;0.023776747,0.049575116,-0.016973166;-0.47405583,0.14382160,-0.063875362;0.30342704,-0.16212738,-0.15749684;-0.010164713,-0.033054799,0.0080676321;-0.068149090,-0.014145531,-0.083623439;0.080996461,0.046955571,0.065057330;-0.80056715,0.61641628,0.27326876;0.0065048584,0.066426113,0.060755216;0.048401710,0.27684277,-0.035768442;-0.37817496,0.46815443,-0.37690166;0.50253326,0.64433789,-0.22256881;0.10804328,-0.10651948,-0.072530910;-0.26477030,0.074723244,0.90906197;0.12373004,-0.10384084,0.10690383;-0.092569664,0.0036151183,-0.060148466;0.067734979,-0.060392909,0.097299635;0.069655597,0.12486014,-0.38743642;-0.047446240,-0.079823971,0.077870943;-0.080286898,0.062817372,0.041379049;-0.11184546,-0.13207248,0.0066199633;-0.10291866,-0.17490597,0.0041896808;-0.61680645,-0.16060983,0.46799311;0.15402365,-0.12539797,-0.23459989;-0.13094640,0.26119620,0.090633251;0.12446557,0.53983581,0.056702100;-0.064821243,0.0028665587,-0.17789589;0.0028397962,0.0022396985,-0.022950187;-0.031183824,-0.034105681,-0.079489440;-0.27047288,0.13264139,0.68019158;-0.18365066,0.30783525,-0.23641501;-0.084395528,0.014986082,0.038748913;0.10412025,0.25519544,-0.21013555;-0.0055523822,-0.18938878,0.13587202;-0.097777925,0.11086731,-0.20151474;0.013454798,0.41120538,-0.29145312;-0.22310972,0.033346836,-0.17319587;0.13977863,-0.078792989,-0.16018809;0.060350846,-0.19806568,0.15751088;0.059629247,-0.12293994,-0.013128475;0.077475309,0.038775988,-0.13875248;-0.073074162,0.034561306,-0.070515536;-0.030371210,0.056076773,-0.0060104541;-0.014176912,-0.023446079,0.096781202;0.17940797,-0.49309134,0.098765910;-0.036261149,0.049778804,0.10537712;-0.22165251,0.26658627,-0.54659730;-0.58360261,-0.27241343,-0.89636153;0.37749493,-0.082415551,0.18595405;0.18782450,-0.21810502,0.29984692;-0.11521876,0.0088482201,-0.25552690;0.15157261,0.091358416,0.044465590;0.027066337,-0.0020436039,0.0099261757;-0.0095810471,-0.12034760,0.042809479;-0.70267987,0.22362821,-0.013574808;-0.00094062253,0.056101121,0.087921351;-0.094699129,0.0052244687,-0.080349661;-0.020348551,-0.036629040,-0.014087967;-0.28562847,0.29433757,0.29889104;0.0085682636,0.031620052,-0.063786954;0.29873917,0.072397731,0.019072488;-0.018752592,-0.73161250,0.26314855;-0.069935955,0.020096350,-0.065628201;-0.93576330,0.10149506,0.091749728;0.11805499,-0.15134864,-0.062733077;-0.036206994,-0.0034468842,0.024932727;-0.19988284,-0.46451011,-0.30865088;0.062685460,-0.0031095541,-0.017761594;0.41865584,-0.29280898,0.24406338;-0.11161800,0.38439572,-0.38964561;0.052199781,0.12316079,-0.038841818;-0.44938046,0.54170841,0.021401649;-0.066964418,0.14784886,0.16460386;0.020565340,0.19741456,0.14665382;-0.15669008,-0.034639638,0.092525035;-0.040721588,0.014005892,-0.076414824;0.14838009,-0.0034469718,0.029192161;-0.011221437,0.026614876,-0.030976078;0.086632401,0.062941872,-0.15563296;-0.00087944325,0.029275060,0.040189557;0.52582347,-0.0029483975,0.079895623;0.078396693,-0.033226706,-0.022581844;-0.068611413,-0.25247151,0.26112196;0.14652574,-0.25601223,0.15691315;-0.094436973,0.43009007,-0.23251858;0.12606563,0.20324616,-0.084253490;-0.14752910,-0.39157853,0.12331896;0.010458344,-0.088903971,0.042526346;0.15649101,-0.026142798,-0.10158617;0.045786846,0.023065135,-0.065181457;0.12724118,-0.10430059,0.032683432;0.077538870,-0.0061167628,0.089893922;0.0088845640,-0.054132938,0.084967442;-0.028500363,0.089753777,0.047505032;0.062406089,-0.35982975,0.051263090;-0.18315062,0.64332592,-0.071338303;0.36296150,-0.77867168,-0.012068376;0.24601854,-0.19535246,0.16144341;-0.080928281,-0.031555668,0.067006573;0.80217355,0.20835005,0.095039114;-0.38920143,0.38995239,-0.092344165;0.066495918,0.067145362,-0.043912802;0.069437757,0.010122486,0.017699068;0.57582808,-0.14772484,-0.22683568;0.44945616,-0.0087714959,-0.49748665;-0.24464531,-0.35656682,0.096572340;0.27354816,0.0091607701,-0.040349800;-0.071957439,0.083961904,0.066977091];
                   obj.name   = 'exp_300';
               case 'comb_300'
                   obj.getDefault();
               otherwise
                   error('Parrameters of this type do not exist.');
           end
       end
       
       function init(obj,varargin)
           % Initializes the random phases
           if isempty(varargin)
               obj.corr_dist = 10;
           else
               obj.corr_dist = varargin{1};
           end
           obj.sos_phase = 2*pi*(rand(obj.Num_sin,length(varargin{1}))-0.5);
       end
       
       function change_corr_dist(obj,corr_dist)
           obj.corr_dist = corr_dist;
           if length(corr_dist) ~= size(obj.sos_phase,2)
               obj.sos_phase = 2*pi*(rand(obj.Num_sin,length(corr_dist))-0.5);
               fprintf('The initial phase was changed.');
           end
       end
       
       function value = randn(obj,pos)
           % input: pos [3 x N]
           pos = repmat(pos,1,1,length(obj.corr_dist)).*repmat(reshape(obj.def_corr_dist./obj.corr_dist,1,1,[]),size(pos,1),size(pos,2),1);
           value = squeeze(multiprod(obj.sos_an,cos(multiprod(obj.sos_wn,pos,[1,2],[1,2]) + permute(repmat(obj.sos_phase,1,1,size(pos,2)),[1,3,2])), [1,2],1));
       end
       
       function value = rand(obj,pos)
           % input: pos [3 x N]
           value = obj.randn(pos);
           % Transform Normal distribution to Uniform distribution
           value = 0.5 * erfc( -value/sqrt(2) );
       end
       
       function value = randi(obj,pos,imax)
           % input: pos [3 x N]
           value = obj.randn(pos);
           % Transform Normal distribution to Uniform distribution
           value = 0.5 * erfc( -value/sqrt(2) );
           value = ceil(value*imax);
       end
       
       %% get and set function
       function out = get.Num_sin(obj)
           out = length(obj.sos_an);
       end
       function out = get.Dim(obj)
           out = size(obj.sos_wn,2);
       end
   end   
end