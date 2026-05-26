(unit=11): JPL00   RATES         (rxns ver 6.0b)   (mjp 3/00)
RATES-----  200  230 RATE COEFFs
RATE/    1        0.               jO2+*=O+O
RATE/    2        0.               jO3+*=O(1D)+O2
RATE/    3        0.               jO3+*=O(ALL)+O2
RATE/1   4  1.80E-11     -110.     O(1D)+M=O+M      --NASA92--
         4  3.30E-11      -55.     O(1D)+O2=O+O2    --JPL06--
RATE/    5  1.63E-10      -60.     O(1D)+H2O=OH+OH  --JPL06--
RATE/    6  1.10E-10        0.     O(1D)+H2=OH+H    --NASA97--
RATE/    7  1.66E-10        0.     O(1D)+CH4=OH+CH3 --JPL09--
RATE/    8  0.09E-10        0.     O1D+CH4=H2+CH2O  --JPL09--
RATE/    9  7.25E-11      -20.     O(1D)+N2O=NO+NO  --JPL09--
RATE/   10  4.63E-11      -20.     O(1D)+N2O=N2+O2  --JPL09--
RATE/   11  6.00E-34       2.4     O+O2+M=O3+M      --NASA00++
RATE/   12  8.00E-12     2060.     O+O3=O2+O2       --NASA92--
RATE/   13  9.59E-34     -480.     O+O+M=O2+M       CAMPBELL+GRAY    XXX
RATE/   14   5.3E-11     5100.     O+H2=OH+H        WESTENBERG       XXX
RATE/   15
RATE/   16  2.20E-11     -120.     O+OH=O2+H        --NASA92--
RATE/   17  3.00E-11     -200.     O+HO2=OH+O2      --NASA00--
RATE/   18  1.40E-12     2000.     O+H2O2=OH+HO2    --NASA92--
RATE/   19  5.10E-12     -210.     O+NO2=NO+O2      --JPL06--
RATE/   20        0.        0.     jNO2+*=NO+O
RATE/   21  3.00E-12     1500.     O3+NO=NO2+O2     --NASA00++
RATE/   22  1.20E-10        0.     O(1D)+O3=O2+O2   --NASA92--
RATE/   23  9.00E-32       1.5     O+NO+M=NO2+M     --NASA92--(k3)
RATE/   24  2.50E-31       1.8     O+NO2+M=NO3+M    --NASA02--(k3)
RATE/   25  1.40E-10      470.     O3+H=OH+O2       --NASA92--
RATE/   26  1.70E-12      940.     O3+OH=HO2+O2     --NASA00++
RATE/   27  1.00E-14      490.     O3+HO2=OH+O2+O2  --NASA02++
RATE/   28  1.20E-13     2450.     O3+NO2=O2+NO3    --NASA92--
RATE/   29  2.10e-11     2200.     O+OCS=CO+SO      --NASA92--
RATE/   30  1.10e-13     1200.     OH+OCS=...         JPL 
RATE/   31        0.        0.     jH2O2+*=OH+OH
RATE/   32  5.70E-32       1.6     H+O2+M=HO2+M     --NASA92--(k3)
RATE/   33  7.20E-11      000.     H+HO2=OH+OH      --JPL06--
RATE/   34  6.90E-12      000.     H+HO2=H2+O2      --JPL06--
RATE/   35  3.50E-37       0.6     O(1D)+N2+M=N2O+M --NASA92--(k3)
RATE/   36  1.00E-11        0.     O+NO3=O2+NO2     --NASA92--
RATE/   37        0.        0.      dens O(1D)
RATE/   38      0.00        0.     <H2O2  source      input:(ppt/d)
RATE/   39  1.80E-12        0.     OH+OH=H2O+O      --JPL06--
RATE/1  40  4.80E-11     -250.     OH+HO2=H2O+O2    --NASA00--
        40        0.        0.       OH+HO2+M  (k2+k3*M)        
RATE/   41  2.90E-12      160.     OH+H2O2=H2O+HO2  --NASA02--
RATE/   42  3.50E-12     -250.     HO2+NO=OH+NO2    --NASA97--
RATE/   43    0.0000        0.     >a:HO2=...       (STICK-PROB)
RATE/2  44  2.30E-13     -600.     HO2+HO2=H2O2+O2  --NASA92--
        44  1.70E-33    -1000.      HO2+HO2+M   (K2 + K3*M)     
        44  2.16E-18        0.      HO2+HO2 AUGMENT BY (1+K*H2O)
RATE/   45       5.0               >r:H2O2/RAIN     (IN DAYS)
RATE/   46    0.0000        0.     >a:H2O2=...      (STICK-PROB)
RATE/   47  2.80E-12     1800.     OH+H2=H2O+H      --NASA92--
RATE/1  48  1.50E-13               OH+CO=CO2+H      --NASA92--
        48  2.44E-20                   =   *(1+.6*P)
RATE/   49  2.45E-12     1775.     OH+CH4=CH3+H2O   --NASA97--
RATE/   50
RATE/   51  3.00E-12     -280.     CH3OO+NO=RO+NO2  --NASA97--
RATE/   52  3.80E-13     -800.     CH3OO+HO2=ROOH+O2--NASA92--
RATE/   53  2.50E-13     -190.     CH3OO+CH3OO=2RO+O2-NASA92--CH3OH=65% PATH?
RATE/   54       0.0               jCH3OOH+*=CH3O+OH
RATE/   55  2.66E-12     -200.     CH3OOH+OH=ROO+H2O -NASA94**70%(R157)
RATE/   56       5.0               >r:CH3OOH/RAIN   (IN DAYS)
RATE/   57    0.0000        0.     >a:CH3OOH=...    (STICK-PROB)
RATE/   58  9.00E-12        0.     H2CO+OH=HCO+H2O  --NASA92--
RATE/   59        0.        0.     jH2CO+*=H+HCO (A)
RATE/   60        0.        0.     jH2CO+*=H2+CO (B)
RATE/   61       5.0               >r:H2CO/RAIN     (IN DAYS)
RATE/2  62  2.00E-30       3.0     NO2+OH=HNO3      --NASA02++(k3)
        62  2.50E-11       0.0      HNO3     (K2)
        62    -.5108        0.      HNO3 EXP=LN(0.6)
RATE/   63        0.        0.     jHNO3+*=OH+NO2
RATE/2  64  2.40E-14     -460.     HNO3+OH=H2O+NO3  --NASA00++(K0)
        64  2.70E-17    -2199.       (K2)
        64  6.50E-34    -1335.       (K3) K=K0+K2*K3*M/(K2+K3*M)
RATE/   65       5.0               >r:HNO3/RAIN=NO     (IN DAYS)
RATE/   66
RATE/2  67  7.00E-31       2.6     NO+OH=HONO       --NASA97--(k3)
        67  3.60E-11       0.1      HONO     (K2)
        67    -.5108        0.      HONO EXP = LN(0.6)
RATE/   68        0.               NO2+HO2=HONO+O2  zeroed out
RATE/   69        0.        0.     jHONO+*=OH+NO
RATE/   70  1.80E-11      390.     HONO+OH=H2O+NO2  --NASA92--
RATE/   71       5.0               >r:HONO/RAIN     (IN DAYS)
RATE/2  72  1.80E-31       3.2     HO2+NO2=HNO4     --NASA92--(k3)
        72  4.70E-12       1.4     HNO4     (K2)
        72    -.5108        0.      HNO4 EXP=LN(0.6)
RATE/   73   4.76E26    10900.     HNO4=HO2+NO2     --NASA92--k(EQ)
RATE/   74        0.        0.     jHNO4+*=OH+NO3
RATE/   75       5.0               >r:HNO4/RAIN      (IN DAYS)
RATE/   76   1.3E-12     -380.     HNO4+OH=H2O+NO2+O2-NASA92--
RATE/   77        0.        0.     jNO3+*=NO2+O
RATE/   78        0.        0.     jNO3+*=NO+O2
RATE/   79  1.50E-11     -170.     NO3+NO=2NO2      --NASA92--
RATE/   80  4.50E-14     1260.     NO3+NO2=NO+O2+NO2--NASA92--
RATE/   81     0.000        0.     <H2CO  source    input:(ppt/d)
RATE/2  82  2.00E-30       4.4     NO2+NO3=N2O5     --NASA00++(k3)
        82  1.40E-12       0.7      N2O5     (K2)
        82    -.5108        0.      N2O5 EXP=LN(0.6)
RATE/   83  3.333E26    10991.     N2O5=NO2+NO3     --NASA00++k(EQ)
RATE/   84        0.        0.     jN2O5+*=NO2+NO3
RATE/   85   0.00E-3        0.     >a:NO2/AEROSOL   (STICK-PROB)
RATE/   86        0.        0.     jNO+*=N+O
RATE/   87  1.50E-11     3600.     N+O2=NO+O        --NASA92--
RATE/   88  2.00E-26               N+O3=NO+O2       --NASA94--LIM/E10
RATE/   89  2.10E-11     -100.     N+NO=N2+O        --NASA94--
RATE/   90  5.80E-12     -220.     N+NO2=N2O+O      --NASA94--
RATE/   91                          dens N(4S)
RATE/   92       1.0 2.00E-20      jC3H6O=2CH3+CO   ---97---/(A+B*M)
RATE/   93  2.20E-12      685.     C3H6O+OH=CH3+CH2O+CO2+H2O --NASA97--
RATE/   94      0.00        0.     <CH3OOH source   input:(ppt/d)
RATE/   95  8.70E-12     1070.     C2H6+OH=...      --NASA92--
RATE/   96      0.00        0.     >r:C2H6    uniform loss
RATE/   97  7.70E-11       90.     C2H6+CL=HCL+...  --NASA92--
RATE/   98
RATE/   99
RATE/  100  2.60E-12      350.     HCL+OH=CL+H2O    --NASA92--
RATE/  101  1.00E-11     3300.     HCL+O=CL+OH      --NASA92--
RATE/  102  1.00E-30        0.     HCL+CLNO3=HNO3+CL2-NASA92--lim/E10
RATE/  103  3.60E-11      375.     CL+HO2=CLO+OH    --JPL09-- (added CAM, 2010)
RATE/  104  9.60E-12     1360.     CL+CH4=HCL+CH3   --NASA00++
RATE/  105        0.        0.     CL+HNO4=HCL+NO2+O2-dropped-
RATE/  106  8.10E-11       30.     CL+H2CO=HCL+HCO  --NASA92--
RATE/  107  3.70E-11     2300.     CL+H2=HCL+H      --NASA92--
RATE/  108  1.40E-11     -269.     CL+HO2=HCL+O2    --JPL09--
RATE/  109  1.10E-11      980.     CL+H2O2=HCL+HO2  --NASA02--
RATE/  110  3.30E-10        0.     O(1D)+CCL4.......--NASA92--
RATE/  111  2.30E-10        0.     O(1D)+CFCL3......--NASA92--
RATE/  112  1.40E-10        0.     O(1D)+CF2CL2=....--NASA92--
RATE/  113  1.50E-10        0.     O(1D)+HCL=OH+CL  --NASA92--
RATE/  114  2.30E-11      200.     CL+O3=CLO+O2     --NASA00++
RATE/  115  2.80E-11      -85.     CLO+O=CL+O2      --JPL06--
RATE/  116  6.40E-12     -290.     CLO+NO=CL+NO2    --NASA92--
RATE/  117        0.        0.     jCLO+*=CL+O
RATE/2 118  1.60E-32       4.5     CLO+CLO=CL2O2    --JPL09++(k3)
       118  3.00E-12       2.0      CL2O2    (K2)
       118    -.5108        0.      CL202EXP=LN(0.6)
RATE/  119  7.874E26     8744.     CL2O2=CLO+CLO    --NASA00++k(EQ)
RATE/  120  6.00E-13     -230.     CLO+OH=HCL+O2    --NASA02++
RATE/  121        0.        0.     CLO+HO2=HCL+O3   --NASA92--(removed CAM, 2010)
RATE/2 122  1.80E-31       3.4     CLO+NO2=CLNO3    --NASA00--(k3)
       122  1.50E-11       1.9      CLNO3    (K2)
       122    -.5108        0.      CLNO3 EXP=LN(0.6)
RATE/  123  2.00E+27    12980.     CLNO3=CLO+NO2    --NASA79--k(EQ) XXX
RATE/  124      0.00        0.     jCLNO3+*=CLO+NO2
RATE/  125      1.00        0.     jCLNO3+*=CL+NO3
RATE/  126  2.90E-12      800.     CLNO3+O=CLO+NO3  --NASA92--
RATE/  127  1.20E-12      330.     CLNO3+OH=HOCL+NO3--NASA92--
RATE/  128  2.60E-12     -290.     CLO+HO2=HOCL+O2  --JPL09--
RATE/  129        0.        0.     jHOCL+*=OH+CL
RATE/  130  3.00E-12      500.     HOCL+OH=H2O+CLO  --NASA92--
RATE/  131                         CH4(OXID)=ODD-H
RATE/  132                         PROD(ODD-H)
RATE/  133                         LOSS(ODD-H)
RATE/  134                         CL-X (loss Ox)
RATE/  135                         PROD(ODD-O)
RATE/  136                         LOSS(ODD-O)
RATE/  137        0.        0.     jN2O+*=N2+O
RATE/  138                         BR-X (loss Ox)
RATE/  139        0.        0.     jCF2CL2+*=...
RATE/  140        0.        0.     jCFCL3+*=...
RATE/  141        0.        0.     jCCL4+*=...
RATE/  142        0.        0.     jCH3CL+*=...
RATE/  143  2.40E-12     1250.     CH3CL+OH=...     --NASA02--
RATE/  144        0.        0.     jCH3CCL3+*=...
RATE/  145  1.80E-12     1550.     CH3CCL3+OH=...   --NASA92--
RATE/  146        0.        0.     jOCS+*=CO+s         
RATE/  147        0.        0.     jCH3BR+*=...
RATE/  148  2.35E-12     1300.     CH3BR+OH=...     --NASA02--
RATE/  149  1.91e-11     -215.     BrNO3+O=BrO+NO3  CMadded,02/06
RATE/  150  1.70e-11     -250.     BrO+OH=Br+HO2    --JPL2006--(CM,11/07)
RATE/  151
RATE/  152
RATE/  153
RATE/  154
RATE/  155
RATE/  156
RATE/  157  1.14E-12     -200.     CH3OOH+OH=CH2O+OH--NASA94**30%(R55)
RATE/  158        0.               jCL2O2+*=CL+CL+O2 
RATE/  159  1.60E-33     -800.     CL+CL+M=CL2+M    --???88--     XXX
RATE/  160        0.               jCL2+*=CL+CL      
RATE/  161   3.5E-13     1370.     CLO+CLO=OCLO+CL  --NASA94--other chan
RATE/  162        0.               jOCLO+*=CLO+O     
RATE/  163  2.50E-12      600.     OCLO+NO=CLO+NO2  --NASA92-- 
RATE/  164  4.50E-13     -800.     OCLO+OH=HOCL+O2  --NASA92-- 
RATE/  165  3.40E-11     -160.     OCLO+CL=CLO+CLO  --NASA92--
RATE/  166  6.50E-12     -135.     CL+CLNO3=CL2+NO3 --NASA97-- 
RATE/  167  2.40E-11        0.     CL+NO3=CLO+NO2   --NASA94** 
RATE/  168     0.000        0.     <O3      source   input:(ppb/d)
RATE/  169        0.        0.     netp:Ox v/v/day
RATE/  170     0.100        0.     >a:N2O5=2*HNO3--  --91-- heterogeneous
RATE/  171    2.E-16    -6200.     >a:CLNO3=HNO3+HOCL-- 0.006exp(-.15(T-200))
RATE/  172  4.00E-37        0.     >a:HCL+CLNO3=CL2+HNO3-- 91--    MJP /91/
RATE/  173  4.00E-37        0.     >a:HCL+N2O5=CL+NO2+HNO3 91--    #170-173
RATE/  174  4.00E-37        0.     >a:HCL+HOCL=CL2+H2O     91--
RATE/  175    0.0000        0.     >a:H2CO=...
RATE/  176    0.0000        0.     >a:NO3=HNO3
RATE/  177    0.0000        0.     >a:NO3=HONO
RATE/  178     0.000        0.     <NO    source      input:(ppt/d)
RATE/  179     0.000        0.     <NO2   source      input:(ppt/d)
RATE/  180        0.        0.     jCFC-113+*=
RATE/  181  2.00E-10        0.     CFC-113+O1D=     --NASA94--
RATE/  182        0.        0.     jCFC-114+*=
RATE/  183  1.30E-10        0.     CFC-114+O1D=     --NASA94--
RATE/  184        0.        0.     CFC-115+*=
RATE/  185  0.50E-10        0.     jCFC-115+O1D=    --NASA94--
RATE/  186        0.        0.     jH-1211+*=    
RATE/  187        0.        0.     jH-1301+*=
RATE/  188        0.        0.     jH-2402+*=
RATE/  189        0.        0.     jHCFC-22+*=
RATE/  190  1.00E-10        0.     HCFC-22+O1D=     --NASA94--
RATE/  191        0.        0.     jHCFC123+*=       
RATE/  192        0.        0.     jHCFC141b+*=
RATE/  193  1.00E-12     1600.     HCFC-22+OH=      --NASA94--
RATE/  194  7.00E-13      900.     HCFC123+OH=      --NASA94--
RATE/  195  1.70E-12     1700.     HCFC141b+OH=     --NASA94--
RATE/  196  1.50E-10        0.     H-1211+O(1D)=    --NASA94--
RATE/  197  1.00E-10        0.     H-1301+O(1D)=    --NASA94--
RATE/  198  1.60E-10        0.     H-2402+O(1D)=    --NASA94--
RATE/  199   3.0E-11     2450.     CLO+CLO=CL+CL+O2 --NASA94--
RATE/  200   1.0E-12     1590.     CLO+CLO=CL2+O2   --NASA94--
RATE/  201  5.50E-12     -200.     HBR+OH=BR+H2O    --JPL06--
RATE/  202  5.80E-12     1500.     HBR+O=BR+OH      --NASA92--
RATE/  203     0.100        0.     >a:BRNO3=HOBr+HNO3--heterogeneous
RATE/  204  4.80E-12      310.     BR+HO2=HBR+O2    --JPL06--
RATE/  205  1.60E-11      780.     BR+O3=BRO+O2     --JPL09---
RATE/  206  1.90E-11     -230.      BRO+O=BR+O2      --NASA97---
RATE/  207  8.80E-12     -260.     BRO+NO=BR+NO2    --NASA92---
RATE/  208  2.00E-27        0.     BRO+O3=BR+2O2    --NASA97--LIM/E10
RATE/1 209  2.40E-12      -40.     BRO+BRO=2BR+O2   --NASA97--
       209  2.80E-14     -860.     BRO+BRO=BR2+O2   --NASA97--
RATE/  210        0.        0.     jBRO+*=BR+O
RATE/  211  4.50E-12     -460.     BRO+HO2=HOBR+O2  --JPL06--
RATE/  212        0.        0.     jHOBR+*=BR+OH
RATE/  213  3.00E-12        0.     HOBR+OH=BRO+H2O  --????85--     XXX
RATE/2 214  5.20E-31       3.2     BRO+NO2=BRNO3    --NASA00++(k3)
       214  6.90E-12       2.9      BRNO3    (K2)
       214    -.5108        0.      BRNO3 EXP=LN(0.6)
RATE/  215       .15        0.     jBRNO3+*=BRO+NO2
RATE/  216       .85        0.     jBRNO3+*=BR+NO3
RATE/  217  4.10E-13     -290.     CLO+BRO=BRCL+O2  --NASA00++
RATE/  218  1.70E-11      800.     BR+H2CO=HBR+HCO  --NASA92--
RATE/  219  9.50E-13     -550.     CLO+BRO=BR+OCLO  --NASA00++
RATE/  220       10.        0.     jBRCL+*=BR+CL    --scale j:CL2
RATE/  221  2.30E-12     -260.     CLO+BRO=BR+CL(O2)--NASA00++
RATE/  222  2.10E-11      830.     I+O3=IO+O2         --JPL19--
RATE/  223  1.40E-10        0.     IO+O=I+O2          --JPL19--
RATE/  224  9.90E-12     -200.     IO+NO=I+NO2        --JPL19--
RATE/  225  1.40E-11     -540.     IO+HO2=HOI+O2      --JPL19--
RATE/2 226  7.70E-31       5.0     IO+NO2+M=IONO2     --JPL19--(k0)
       226  1.60E-11       0.0     IONO2              (kinf)
       226    -.5108        0.     IONO2 EXP=LN(0.6)
RATE/  227  1.50E-11     1090.     I+HO2=HI+O2        --JPL19--
RATE/  228  3.00E-11     1120.     HI+OH=I+H2O        --JPL19--
RATE/  229  5.10E-12     -280.     IO+CLO=CL+I+O2     --JPL19--
RATE/  230  1.50E-11        0.     IO+BRO=BR+I+O2     --JPL19--
