function [ensoFit,ensoRes,ensoData]=ensoFitRes(b)

D=[12.90000    1.000000
    11.30000    2.000000
    10.60000    3.000000
    11.20000    4.000000
    10.90000    5.000000
    7.500000    6.000000
    7.700000    7.000000
    11.70000    8.000000
    12.90000    9.000000
    14.30000   10.000000
    10.90000    11.00000
    13.70000    12.00000
    17.10000    13.00000
    14.00000    14.00000
    15.30000    15.00000
    8.500000    16.00000
    5.700000    17.00000
    5.500000    18.00000
    7.600000    19.00000
    8.600000    20.00000
    7.300000    21.00000
    7.600000    22.00000
    12.70000    23.00000
    11.00000    24.00000
    12.70000    25.00000
    12.90000    26.00000
    13.00000    27.00000
    10.90000    28.00000
   10.400000    29.00000
   10.200000    30.00000
    8.000000    31.00000
    10.90000    32.00000
    13.60000    33.00000
   10.500000    34.00000
    9.200000    35.00000
    12.40000    36.00000
    12.70000    37.00000
    13.30000    38.00000
   10.100000    39.00000
    7.800000    40.00000
    4.800000    41.00000
    3.000000    42.00000
    2.500000    43.00000
    6.300000    44.00000
    9.700000    45.00000
    11.60000    46.00000
    8.600000    47.00000
    12.40000    48.00000
   10.500000    49.00000
    13.30000    50.00000
   10.400000    51.00000
    8.100000    52.00000
    3.700000    53.00000
    10.70000    54.00000
    5.100000    55.00000
   10.400000    56.00000
    10.90000    57.00000
    11.70000    58.00000
    11.40000    59.00000
    13.70000    60.00000
    14.10000    61.00000
    14.00000    62.00000
    12.50000    63.00000
    6.300000    64.00000
    9.600000    65.00000
    11.70000    66.00000
    5.000000    67.00000
    10.80000    68.00000
    12.70000    69.00000
    10.80000    70.00000
    11.80000    71.00000
    12.60000    72.00000
    15.70000    73.00000
    12.60000    74.00000
    14.80000    75.00000
    7.800000    76.00000
    7.100000    77.00000
    11.20000    78.00000
    8.100000    79.00000
    6.400000    80.00000
    5.200000    81.00000
    12.00000    82.00000
   10.200000    83.00000
    12.70000    84.00000
   10.200000    85.00000
    14.70000    86.00000
    12.20000    87.00000
    7.100000    88.00000
    5.700000    89.00000
    6.700000    90.00000
    3.900000    91.00000
    8.500000    92.00000
    8.300000    93.00000
    10.80000    94.00000
    16.70000    95.00000
    12.60000    96.00000
    12.50000    97.00000
    12.50000    98.00000
    9.800000    99.00000
    7.200000   100.00000
    4.100000   101.00000
    10.60000   102.00000
   10.100000   103.00000
   10.100000   104.00000
    11.90000   105.00000
    13.60000    106.0000
    16.30000    107.0000
    17.60000    108.0000
    15.50000    109.0000
    16.00000    110.0000
    15.20000    111.0000
    11.20000    112.0000
    14.30000    113.0000
    14.50000    114.0000
    8.500000    115.0000
    12.00000    116.0000
    12.70000    117.0000
    11.30000    118.0000
    14.50000    119.0000
    15.10000    120.0000
   10.400000    121.0000
    11.50000    122.0000
    13.40000    123.0000
    7.500000    124.0000
   0.6000000    125.0000
   0.3000000    126.0000
    5.500000    127.0000
    5.000000    128.0000
    4.600000    129.0000
    8.200000    130.0000
    9.900000    131.0000
    9.200000    132.0000
    12.50000    133.0000
    10.90000    134.0000
    9.900000    135.0000
    8.900000    136.0000
    7.600000    137.0000
    9.500000    138.0000
    8.400000    139.0000
    10.70000    140.0000
    13.60000    141.0000
    13.70000    142.0000
    13.70000    143.0000
    16.50000    144.0000
    16.80000    145.0000
    17.10000    146.0000
    15.40000    147.0000
    9.500000    148.0000
    6.100000    149.0000
   10.100000    150.0000
    9.300000    151.0000
    5.300000    152.0000
    11.20000    153.0000
    16.60000    154.0000
    15.60000    155.0000
    12.00000    156.0000
    11.50000    157.0000
    8.600000    158.0000
    13.80000    159.0000
    8.700000    160.0000
    8.600000    161.0000
    8.600000    162.0000
    8.700000    163.0000
    12.80000    164.0000
    13.20000    165.0000
    14.00000    166.0000
    13.40000    167.0000
    14.80000    168.0000];

y=D(:,1);
x=D(:,2);

res=b(1) + b(2)*cos( 2*pi*x/12 ) + b(3)*sin( 2*pi*x/12 ) + b(5)*cos( 2*pi*x*b(4) ) ;
res=res+b(6)*sin( 2*pi*x*b(4) ) + b(8)*cos( 2*pi*x*b(7) ) + b(9)*sin( 2*pi*x*b(7) );
ensoFit=res;
ensoRes=y-res;
ensoData=y;

