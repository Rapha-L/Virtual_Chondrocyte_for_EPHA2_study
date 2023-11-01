function []=Save_pie_chart(condition)

%% This file is part of the repository Virtual_Chondrocyte_EPHA2_study
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
%%

%Titles= ["EPHA2 up","Inflammation", "Inflammation + EPHA2 up", "Inflammation + EPHA2 down"];
p_labels = ["Healthy","Hypertrophic"];
title_size =11;

%     p = pie(condition(c).result_transition);
%     legend(labels)
%     saveas(p,fullfile('Figures',['Infl_nodes.fig'))
%     
t = tiledlayout(2,3);

ax1 = nexttile;
pie(ax1,[100,0])
lgd=legend(p_labels);
lgd.FontSize = 12;
title({"Initial state";''},'FontSize',title_size)

c=1;
ax2 = nexttile;
pie(ax2,condition(c).result_transition)
title({"EPHA2 +";''},'FontSize',title_size)

c=2;
ax3 = nexttile;
pie(ax3,condition(c).result_transition)
title({"Inflammation +";''},'FontSize',title_size)

c=3;
ax4 = nexttile;
pie(ax4,condition(c).result_transition)
title({'Inflammation +';' EPHA2 +';''},'FontSize',title_size)

c=4;
ax5 = nexttile;
pie(ax5,condition(c).result_transition)
title({'Inflammation +';' EPHA2 -';''},'FontSize',title_size)

end