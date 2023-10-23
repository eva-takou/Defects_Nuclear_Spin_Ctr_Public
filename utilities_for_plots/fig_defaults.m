function fig_defaults(FntSize)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Set default properties for plots.
%--------------------------------------------------------------------------

set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.01 .01], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1)

set(gcf,'color','w')
set(gca,'fontsize',FntSize,'fontname','Microsoft Sans Serif')
grid off
end