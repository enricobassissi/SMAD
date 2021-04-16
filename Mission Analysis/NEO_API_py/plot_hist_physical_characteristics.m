function h_fig_pha_char_dist = plot_hist_physical_characteristics(physical_properties,colors,ps_vector)

    h_fig_pha_char_dist = figure('Name','PHAs Characteristics Distribution');
    perc_labels = {'25 Perc','50 Perc','75 Perc'};%{'Sample','x=12'}
    % title('PHAs Characteristics Distribution')
    subplot(2,2,1)
    % h1=histogram(physical_properties(:,1));
    h1=histogram(physical_properties(:,1));%weibull
    h1.FaceColor=colors(1,:);
    h1.EdgeColor=colors(1,:);
    h1.FaceAlpha=1;
    morebins(h1); morebins(h1); morebins(h1); morebins(h1);
    % h1(2).Color=colors(2,:);
    % yscale logarithm
    xlabel('(a)')
    ylabel('Asteroid Count (log)')
    legend('MOID (AU)')
    % set(gca,'XScale','log')
    set(gca,'YScale','log')
    Y_1 = prctile(physical_properties(:,1),[25,50,75],1);
    hold on
    xl(1,1) = xline(Y_1(1),'-','Color',colors(4,:),'LineWidth', 1.6,'DisplayName',perc_labels{1});
    xl(1,2) = xline(Y_1(2),'-.','Color',colors(4,:),'LineWidth', 1.6,'DisplayName',perc_labels{2});
    xl(1,3) = xline(Y_1(3),'--','Color',colors(4,:),'LineWidth', 1.6,'DisplayName',perc_labels{3});
    % xl4(1, 1).Label = '25 Perc';
    legend('show')
    legend('Location','northeast');

    subplot(2,2,2)
    histogram(physical_properties(:,2),'FaceColor',colors(1,:),'FaceAlpha',1);
    xlabel('(b)')
    ylabel('Asteroid Count ($\#$)')
    legend('OCC')
    Y_2 = prctile(physical_properties(:,2),[25,50,75],1);
    hold on
    xl(2,1) = xline(Y_2(1),'-','Color',colors(4,:),'LineWidth', 1.6,'DisplayName',perc_labels{1});
    xl(2,2) = xline(Y_2(2),'-.','Color',colors(4,:),'LineWidth', 1.6,'DisplayName',perc_labels{2});
    xl(2,3) = xline(Y_2(3),'--','Color',colors(4,:),'LineWidth', 1.6,'DisplayName',perc_labels{3});
    % xl4(1, 1).Label = '25 Perc';
    legend('show')
    legend('Location','northwest');

    subplot(2,2,3)
%     h3=histfit(physical_properties(:,3),40,'kernel'); % weibull
    h3=histogram(physical_properties(:,3)); % weibull
    h3.FaceColor=colors(1,:);
    h3.EdgeColor=colors(1,:);
    h3.FaceAlpha = 1;
    morebins(h3); morebins(h3); morebins(h3);
%     h3(2).Color=colors(2,:);
    xlabel('(c)')
    ylabel('Asteroid Count ($\#$)')
    legend('H')
    Y_3 = prctile(physical_properties(:,3),[25,50,75],1);
    hold on
    xl(3,1) = xline(Y_3(1),'-','Color',colors(4,:),'LineWidth', 1.6,'DisplayName',perc_labels{1});
    xl(3,2) = xline(Y_3(2),'-.','Color',colors(4,:),'LineWidth', 1.6,'DisplayName',perc_labels{2});
    xl(3,3) = xline(Y_3(3),'--','Color',colors(4,:),'LineWidth', 1.6,'DisplayName',perc_labels{3});
    % xl4(1, 1).Label = '25 Perc';
    legend('show')
    legend('Location','northwest');

    subplot(2,2,4)
%     h4=histfit(ps_vector,400,'kernel');%,'DisplayName','PS'
    h4=histogram(ps_vector);%,'DisplayName','PS'
    legend('PS') %,'Interpolation'
    xlim([-13,0])
    h4.FaceColor=colors(1,:);
    h4.FaceAlpha = 1;
    h4.EdgeColor=colors(1,:);
    morebins(h4); morebins(h4); morebins(h4); morebins(h4); morebins(h4); morebins(h4); morebins(h4);
    morebins(h4); morebins(h4); morebins(h4); morebins(h4); morebins(h4); morebins(h4); morebins(h4);
%     h4(2).Color=colors(2,:);
    xlabel('(d)')
    ylabel('Asteroid Count ($\#$)')
    Y_4 = prctile(ps_vector,[25,50,75],2);
    hold on
    xl(4,1) = xline(Y_4(1),'-','Color',colors(4,:),'LineWidth', 1.6,'DisplayName',perc_labels{1});
    xl(4,2) = xline(Y_4(2),'-.','Color',colors(4,:),'LineWidth', 1.6,'DisplayName',perc_labels{2});
    xl(4,3) = xline(Y_4(3),'--','Color',colors(4,:),'LineWidth', 1.6,'DisplayName',perc_labels{3});
    % xl4(1, 1).Label = '25 Perc';
    legend('show');

end