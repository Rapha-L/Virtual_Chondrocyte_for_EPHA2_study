function [ y ] = zupdates_AC_EPHA2(z)
%% Virtual chondrocyte model (equations) with EPHA2 for SYNCHRONOUS time steps

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
% Ras IGF-IR Msx2 EF-1 ATF4 HIF-2alpha GREM1 FRZB DKK1 Frizzled-LRP5/7
% 41   42    43   44   45     46         47   48   49      50
% Cytokines  ALK1 ALK5 Rinflam TAK1 JNK Proteoglycans IkB-a Socs FOXO1
%   51        52   53    54     55   56     57          58    59   60
% IGFBP EPHA2
%   61   62
%%

global s
global n
y = zeros(3,n);

% switch node
%     case 1
        y(2,1) = 1 - z(1,49) - z(1,48); 
%     case 2
        y(2,2) = z(1,50) + z(1,62) - 0.3*z(1,4); 
%     case 3
        y(2,3) = 1 - z(1,61);
%     case 4
        y(2,4) = z(1,52) - 0.5*z(1,25) + z(1,52); 
%     case 5
        y(2,5) = 1; 
%     case 6
        y(2,6) =  (2/3)*z(1,5);
%     case 7
        y(2,7) = 1 - z(1,36);
%    case 8
        y(2,8) = z(1,7)- z(1,10);
%     case 9
        y(2,9) = z(1,38) + z(1,22) + z(1,32) - z(1,10) - z(1,14)- z(1,26)*z(1,30) - z(1,31) - z(1,43);
%     case 10
        y(2,10) = z(1,14) + z(1,26) - z(1,7) - z(1,9);
%     case 11
        y(2,11) = 1;
%     case 12
        y(2,12) = z(1,11);
%     case 13
        y(2,13) = 1; 
%     case 14
        y(2,14) = (z(1,12) + z(1,3))*(s(2)); 
%     case 15
        y(2,15) = z(1,34) + z(1,19) - z(1,30) - z(1,26);
%     case 16
        y(2,16) = 1; 
%     case 17
        y(2,17) = z(1,16);
%     case 18
        y(2,18) = z(1,17) + z(1,42) - z(1,19) + z(1,54)*(1-z(1,17)) -z(1,59); %
%     case 19
        y(2,19) = z(1,4) - 0.25 * z(1,22);  
%     case 20
        y(2,20) = 1;
%     case 21
        y(2,21) = 1;
%     case 22
        y(2,22) = s(3) *(z(1,41) - z(1,37) + z(1,55) + z(1,14) + z(1,62));
%     case 23
        y(2,23) = 1;
%     case 24
        y(2,24) = 1; 
%     case 25
        y(2,25) = 1;
%     case 26
        y(2,26) = z(1,53) - 0.5*(s(3)*(z(1,22)+ z(1,31)+ z(1,25)));
%     case 27
        y(2,27) = z(1,16);
%     case 28
        y(2,28) = (z(1,34)*z(1,23)*z(1,55) + z(1,34)* z(1,26) + z(1,56))*s(3); 
%     case 29
        y(2,29) = (z(1,27) + z(1,38) + z(1,55) - z(1,58))*s(2); 
%     case 30
        y(2,30) = z(1,37) - z(1,35);
%     case 31
        y(2,31) = 1 - 0.5*(z(1,35) + z(1,18));
%     case 32
        y(2,32) = z(1,34) + z(1,19) - z(1,43);
%     case 33
        y(2,33) = 1-z(1,47);
%     case 34
        y(2,34) = (z(1,23) + z(1,33) + z(1,55))*s(3);
%     case 35
        y(2,35) = 1 - z(1,14) - z(1,38);
%     case 36
        y(2,36) = (1 + z(1,37) - z(1,2) - 0.5*z(1,22));    
%     case 37
        y(2,37) = z(1,14);
%     case 38
        y(2,38) = z(1,39) - 0.5*z(1,37);
%     case 39
        y(2,39) = (z(1,41) + z(1,42) + z(1,54) + z(1,62))*s(4); 
%     case 40
        y(2,40) = s(2) *(z(1,22) + z(1,56));
%     case 41
        y(2,41) = (z(1,1) + z(1,33) +  z(1,17) + z(1,27))*s(4); 
%     case 42
        y(2,42) = z(1,3) - 0.5*z(1,59);
%     case 43
        y(2,43) = 1 - z(1,32);
%     case 44
        y(2,44) = 1;
%     case 45
        y(2,45) = (z(1,14) + z(1,22))*s(2);
%     case 46
        y(2,46) = 1;
%    case 47        
        y(2,47) = 1;
%    case 48        
        y(2,48)= 1;
%    case 49        
        y(2,49)= 1;
%    case 50       
        y(2,50)= z(1,1) - z(1,49);
%   case 51       
        y(2,51)= 1;
%    case 52       
        y(2,52)= z(1,23) + z(1,33);
%    case 53       
        y(2,53)= z(1,23);                
%    case 54       
        y(2,54)= z(1,51) - z(1,59);
%    case 55       
        y(2,55)= z(1,54) - 0.5*z(1,34) + z(1,52);
%    case 56       
        y(2,56)= z(1,55);     
%    case 57       
        y(2,57)= 1 - z(1,24);   
%    case 58       
        y(2,58)= 1 - z(1,55);
%    case 59       
        y(2,59)= 1;        
%    case 60       
        y(2,60)= 0.75 + z(1,56) - z(1,38) - z(1,22) ;
        
        y(2,61) = 1;
        
        y(2,62) = s(2) *(z(1,38) + z(1,1)); 
        
% solw interactions
%     case 47
        y(3,1) = 2*z(1,6) - z(1,43);
%     case 48
        y(3,2) = 1;
%     case 49
        y(3,3) = 1; 
%     case 50
        y(3,4) = 1 - z(1,45) + z(1,44);
%     case 51
         y(3,5) = (z(1,9) + z(1,19) + z(1,29) + z(1,45)  - z(1,44) - z(1,17) -z(1,26) + z(1,43))*s(2); 
%     case 52
        y(3,6) = 1 - z(1,17); 
%     case 53
        y(3,7) = 1;
%     case 54
        y(3,8) = (1 + z(1,9)+ z(1,7))*s(3);
%     case 55
        y(3,9) = (z(1,46) + z(1,8) + z(1,9) + z(1,15) + z(1,32) - z(1,21)*z(1,19) - z(1,43) + z(1,6) + z(1,40))*s(4);
%     case 56
        y(3,10) = (z(1,14)*z(1,45) - z(1,29) + z(1,34)+ z(1,29) + z(1,10) + z(1,21))*s(4);
%     case 57
        y(3,11) = (z(1,6)*z(1,10) + z(1,10)+ z(1,26))*s(3);
%     case 58
        y(3,12) = (z(1,6) + z(1,10) + z(1,19) - z(1,18))*s(2);
%     case 59
        y(3,13) = (z(1,9) + z(1,4) + z(1,15) - z(1,14) - z(1,11) + z(1,46) - z(1,18) + z(1,7))*s(2);
%     case 60
        y(3,14) = 1;
%     case 61
        y(3,15) = (z(1,9) + z(1,19))*s(2);
%     case 62
        y(3,16) = (z(1,7) + z(1,9))*s(2);
%     case 63
        y(3,17) = z(1,10) - z(1,40)*z(1,22); 
%     case 64
        y(3,18) = 1;
%     case 65
        y(3,19) = 1;
%     case 66
        y(3,20) = (z(1,10) + z(1,26) + z(1,21) - z(1,7) + z(1,60))*s(3);
%     case 67
        y(3,21) = (z(1,10) + z(1,14)*z(1,45))*s(2); 
%     case 68
        y(3,22) = z(1,15);
%     case 69
        y(3,23) = z(1,6) - z(1,40) + z(1,26);
%     case 70
        y(3,24) = (z(1,9) + z(1,29) +z(1,46) + z(1,6) + z(1,46)*z(1,34) + z(1,19) + z(1,40) - z(1,26) + z(1,40)*z(1,28)*z(1,56) - z(1,60))/5;% also activated by HIF-2alpha (hirata, 2012) % changed /3 by /(3/2) % changed expression to be maximally sensitive, added C/EBPb: runx2 and p38 activation, hif transcription 
%     case 71
        y(3,25) = (z(1,18) + z(1,19) + z(1,26) + 2*z(1,29))*s(5);
%     case 72
        y(3,26) = 1;
%     case 73
        y(3,27) = z(1,9) - 0.5*z(1,19) +z(1,29);
%     case 74
        y(3,28) = 1; 
%     case 75
        y(3,29) = 1;
%     case 76
        y(3,30) = 1;
%     case 77
        y(3,31) = (z(1,28) + z(1,6) + 1.5*z(1,14) + z(1,60)*z(1,26))*s(4);
%     case 78
        y(3,32) = (z(1,15) + z(1,19) - z(1,26));
%     case 79
        y(3,33) = (z(1,6) + z(1,29))*s(2);
%     case 80
        y(3,34) = 1 - 0.4*z(1,44);
%     case 81
        y(3,35) = 1;
%     case 82
        y(3,36) = 1;
%     case 83
        y(3,37) = (1 - z(1,7))*s(2);
%     case 84
        y(3,38) = z(1,9);
%     case 85
        y(3,39) = z(1,9);
%     case 86
        y(3,40) = 1;
%     case 87
        y(3,41) = 1;
%     case 88
        y(3,42)= (z(1,26) + z(1,29))*s(2);
%     case 89
        y(3,43) = (z(1,19) + z(1,26) + z(1,43) - z(1,32))*s(2);
%     case 90
        y(3,44) = (z(1,29) + z(1,40) - z(1,19) + z(1,4))*s(2);
%     case 91
        y(3,45) = 1;
%     case 92
        y(3,46) = z(1,29); 
%    case 97      
        y(3,47) = z(1,10) - z(1,7);
%    case 98       
        y(3,48) = 0.5*z(1,19) + z(1,45) - z(1,7) + z(1,10) -z(1,28);
%    case 99        
        y(3,49) = 0.5*z(1,19) - z(1,7) + z(1,10) -z(1,29) -z(1,28); 
%    case 100        
        y(3,50) = 1;
%    case n+51        
        y(3,51) = z(1,40)*z(1,28)*z(1,56) + z(1,29);        
%    case n+52        
        y(3,52) = s(3)*(0.4 + z(1,8) - z(1,26) + z(1,9) + z(1,28)*z(1,56));  
%    case n+53        
        y(3,53) = s(2)*(0.4 + z(1,26) - z(1,29) + z(1,10));
%    case n+54        
        y(3,54) = 1 - s(3) *(z(1,60) + z(1,26));        
%    case n+55        
        y(3,55) = 1;        
%    case n+56        
        y(3,56) = 1;        
%    case n+57        
        y(3,57) = s(3)*(z(1,10) + z(1,60) - z(1,29));        
%    case n+58        
        y(3,58) = z(1,29);    
%    case n+59        
        y(3,59) = s(2) *(z(1,42) + z(1,17));
%    case n+60         
        y(3,60) = 0.75 + z(1,23)*z(1,34)*z(1,55) - z(1,40)*z(1,28) - z(1,29) - z(1,56)*z(1,40)*z(1,28) ;
        
        y(3,61) = z(1,9);
        
        y(3,62) = z(1,29) + z(1,41);
  % end

y(y > 1) = 1;
y(y < 0) = 0;
y(1,:) = y(2,:).*y(3,:);

end

