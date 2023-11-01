function [ y ] = zupdateas_AC_EPHA2( z,node,speedstring )
%% Virtual chondrocyte model (equations) with EPHA2 for ASYNCHRONOUS time steps

% This file is part of the repository Virtual_Chondrocyte_EPHA2_study
% (https://github.com/Rapha-L) 

%Copyright (c) 2017-2023 - KU Leuven

%File author(s): R. Lesage

%Distributed under the GPLv3 License.
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>


%% Each gene/protein is associated with an index number:

% Wnt Dsh IGFI  R-Smad Ihh Gli2 Bcat Lef/Tcf Runx2 Sox9 PTHrP PPR Col-X PKA
%  1   2    3     4     5   6    7     8       9    10    11   12  13    14
% MEF2C FGF FGFR3 STAT1 Smadcomplex Col-II Nkx3.2 ERK1/2 TGFbeta 
%   15   16   17    18       19       20     21     22      23
% MMP13 Smad7 Smad3 FGFR1 ATF2 NFkb HDAC4 CCND1 Dlx5  BMP  p38
%   24   25    26    27    28   29   30    31    32   33    34
%   GSK DC PP2A Akt PI3K Elk-1 
%   35  36  37  38   39   40   
% Ras IGF-IR Msx2 EF-1 ATF4 HIF-2alpha  GREM1 FRZB DKK1 Frizzled-LRP5/7
% 41   42    43   44   45     46         47   48   49      50     
% Cytokines  ALK1 ALK5 Rinflam TAK1 JNK Proteoglycans IkB-a Socs FOXO1
%   51        52   53    54     55   56     57          58   59    60
% IGFBP EPHA2
%   61   62
%%

global s
global n
if (strcmp(speedstring,'slow'))
node = node + n;
end
switch node
    case 1
        y = 1 - z(1,49) - z(1,48); 
    case 2
        y = z(1,50) + z(1,62) - 0.3*z(1,4);
    case 3
        y = 1 - z(1,61);
    case 4
        y = z(1,52) - 0.5*z(1,25) + z(1,52);  
    case 5
        y = 1; 
    case 6
        y =  (2/3)*z(1,5); 
    case 7
        y = 1 - z(1,36);
    case 8
        y = z(1,7)- z(1,10);
    case 9
        y = z(1,38) + z(1,22) + z(1,32) - z(1,10) - z(1,14)- z(1,26)*z(1,30) - z(1,31) - z(1,43);
    case 10
        y = z(1,14) + z(1,26) - z(1,7) - z(1,9); 
    case 11
        y = 1; 
    case 12
        y = z(1,11);
    case 13
        y = 1; 
    case 14
        y = (z(1,12) + z(1,3))*(s(2));
    case 15
        y = z(1,34) + z(1,19) - z(1,30) - z(1,26);
    case 16
        y = 1; 
    case 17
        y = z(1,16);
    case 18
        y = z(1,17) + z(1,42) - z(1,19) + z(1,54)*(1-z(1,17)) -z(1,59) ;
        y = z(1,4) - 0.25 * z(1,22);  
    case 20
        y = 1;                                                                         
    case 21
        y = 1;
    case 22
        y = s(3) *(z(1,41) - z(1,37) + z(1,55) + z(1,14) + z(1,62));
    case 23
        y = 1;
    case 24
        y = 1; 
    case 25
        y = 1;
    case 26
        y = z(1,53) - 0.5*(s(3)*(z(1,22)+ z(1,31)+ z(1,25)));
    case 27
        y = z(1,16);   
    case 28
        y = (z(1,34)*z(1,23)*z(1,55) + z(1,34)* z(1,26) + z(1,56))*s(3);
    case 29
        y = (z(1,27) + z(1,38) + z(1,55) - z(1,58))*s(2);
    case 30
        y = z(1,37) - z(1,35);
    case 31
        y = 1 - 0.5*(z(1,35) + z(1,18));
    case 32
        y = z(1,34) + z(1,19) - z(1,43); 
    case 33
        y = 1-z(1,47);
    case 34
        y = (z(1,23) + z(1,33) + z(1,55))*s(3); 
    case 35 %free GSK3B
        y = 1 - z(1,14) - z(1,38);        
    case 36 % GSK3B insulated in DC
        y = (1 + z(1,37) - z(1,2) - 0.5*z(1,22));      
    case 37
        y = z(1,14);
    case 38
        y = z(1,39) - 0.5*z(1,37);
    case 39
        y = (z(1,41) + z(1,42) + z(1,54) + z(1,62))*s(4); 
    case 40
        y = s(2) *(z(1,22) + z(1,56));
    case 41
        y = (z(1,1) + z(1,33) +  z(1,17) + z(1,27))*s(4);
    case 42
        y = z(1,3) - 0.5*z(1,59);
    case 43
        y = 1 - z(1,32);
    case 44
        y = 1;
    case 45
        y = (z(1,14) + z(1,22))*s(2);
    case 46
        y = 1;
    case 47        
        y = 1;
    case 48        
        y= 1;
    case 49        
        y= 1;
    case 50       
        y= z(1,1) - z(1,49);       
    case 51       
        y= 1;
    case 52       
        y= z(1,23) + z(1,33);
    case 53       
        y= z(1,23);                
    case 54       
        y= z(1,51) - z(1,59);
    case 55       
        y= z(1,54) - 0.5*z(1,34) + z(1,52);    
    case 56       
        y= z(1,55);     
    case 57       
        y= 1 - z(1,24);     
    case 58       
        y= 1 - z(1,55);
    case 59
        y= 1;
    case 60
        y= 0.75 + z(1,56) - z(1,38) - z(1,22) ;
    case 61
        y= 1;
    case 62
        y = s(2) *(z(1,38) + z(1,1)); 
        
  % slow variables (slow interactions)     
    case n+1
        y = 2*z(1,6) - z(1,43);
    case n+2
        y = 1;
    case n+3
        y = 1; 
    case n+4
        y = 1 - z(1,45) + z(1,44);
    case n+5
         y = (z(1,9) + z(1,19) + z(1,29) + z(1,45)  - z(1,44) - z(1,17) -z(1,26) + z(1,43))*s(2);      
    case n+6
        y = 1 - z(1,17); 
    case n+7
        y = 1;
    case n+8
        y = (1 + z(1,9)+ z(1,7))*s(3); 
    case n+9
        y = (z(1,46) + z(1,8) + z(1,9) + z(1,15) + z(1,32) - z(1,21)*z(1,19) - z(1,43) + z(1,6) + z(1,40))*s(4); 
    case n+10
        y = (z(1,14)*z(1,45) - z(1,29) + z(1,34)+ z(1,29) + z(1,10) + z(1,21))*s(4);
    case n+11
        y = (z(1,6)*z(1,10) + z(1,10)+ z(1,26))*s(3); 
    case n+12
        y = (z(1,6) + z(1,10) + z(1,19) - z(1,18))*s(2);
    case n+13
        y = (z(1,9) + z(1,4) + z(1,15) - z(1,14) - z(1,11) + z(1,46) - z(1,18) + z(1,7))*s(2); 
    case n+14
        y = 1;
    case n+15
        y = (z(1,9) + z(1,19))*s(2); 
    case n+16
        y = (z(1,7) + z(1,9))*s(2); 
    case n+17
        y = z(1,10) - z(1,40)*z(1,22); 
    case n+18
        y = 1;
    case n+19
        y = 1;
    case n+20
        y = (z(1,10) + z(1,26) + z(1,21) - z(1,7) + z(1,60))*s(3); 
    case n+21
        y = (z(1,10) + z(1,14)*z(1,45))*s(2); 
    case n+22
        y = z(1,15);
    case n+23
        y = z(1,6) - z(1,40) + z(1,26);
    case n+24
        y = (z(1,9) + z(1,29) +z(1,46) + z(1,6) + z(1,46)*z(1,34) + z(1,19) + z(1,40) - z(1,26) + z(1,40)*z(1,28)*z(1,56) - z(1,60))/5;
    case n+25
        y = (z(1,18) + z(1,19) + z(1,26) + 2*z(1,29))*s(5); 
    case n+26
        y = 1;
    case n+27
        y = z(1,9) - 0.5*z(1,19) +z(1,29);
    case n+28
        y = 1;
    case n+29
        y = 1;
    case n+30
        y = 1;
    case n+31
        y = (z(1,28) + z(1,6) + 1.5*z(1,14) + z(1,60)*z(1,26))*s(4); 
    case n+32
        y = (z(1,15) + z(1,19) - z(1,26));
    case n+33
        y = (z(1,6) + z(1,29))*s(2); 
    case n+34
        y = 1 - 0.4*z(1,44);
    case n+35
        y = 1;
    case n+36
        y = 1;
    case n+37
        y = (1 - z(1,7))*s(2);
    case n+38
        y = z(1,9);
    case n+39
        y = z(1,9);
    case n+40
        y = 1;
    case n+41
        y = 1;
    case n+42
        y = (z(1,26) + z(1,29))*s(2);
    case n+43
        y = (z(1,19) + z(1,26) + z(1,43) - z(1,32))*s(2); 
    case n+44
        y = (z(1,29) + z(1,40) - z(1,19) + z(1,4))*s(2); 
    case n+45
        y = 1;
    case n+46
        y = z(1,29);       
     case n+47      
        y = z(1,10) - z(1,7);
    case n+48       
        y = 0.5*z(1,19) + z(1,45) - z(1,7) + z(1,10) -z(1,28);
    case n+49        
        y = 0.5*z(1,19) - z(1,7) + z(1,10) -z(1,29) -z(1,28);
    case n+50        
        y = 1;
    case n+51        
        y = z(1,40)*z(1,28)*z(1,56) + z(1,29);        
    case n+52        
        y = s(3)*(0.4 + z(1,8) - z(1,26) + z(1,9) + z(1,28)*z(1,56));
    case n+53        
        y = s(2)*(0.4 + z(1,26) - z(1,29) + z(1,10));
    case n+54        
        y = 1 - s(3) *(z(1,60) + z(1,26));        
    case n+55        
        y = 1;        
    case n+56        
        y = 1;        
    case n+57        
        y = s(3)*(z(1,10) + z(1,60) - z(1,29));        
    case n+58        
        y = z(1,29);    
    case n+59        
        y = s(2) *(z(1,42) + z(1,17));
    case n + 60
        y = 0.75 + z(1,23)*z(1,34)*z(1,55) - z(1,40)*z(1,28) - z(1,29) - z(1,56)*z(1,40)*z(1,28) ;
    case n+ 61
        y = z(1,9);
    case n + 62
        y = z(1,29) + z(1,41);
end

y(y > 1) = 1; 
y(y < 0) = 0;


end

