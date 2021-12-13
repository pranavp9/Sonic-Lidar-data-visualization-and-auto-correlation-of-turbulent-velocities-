clear all
clc

load ('MAE_577_sonic_segment.mat')

k=1;

% U mean for 15 min
for i=0:15*600:252000-(15*600)
    avg_time_1_min(k,1)=mean(time([i+1:i+15*600,1]));
    k=k+1;
end

% k value change so re-intialize k=1; or change k variable to something
% else
k=1;

% time mean for 15 min
for i=0:600*15:252000-(15*600)
    avg_u_1_min(k,1)=mean(u([i+1:i+15*600,1]));
    k=k+1;
end

% u bar
ubar = mean(avg_u_1_min);
ubar_vec_array = zeros(length(avg_time_1_min),1);
ubar_vec_array(:)= ubar;

% plot of u bar over 15 min
hold on
figure(1)
plot(avg_time_1_min,avg_u_1_min)
% mean u bar over time
%plot(avg_time_1_min,ubar_vec_array)
title("u bar Vs time running mean 15 min")
xlabel('Average time 15 min')
ylabel("bar(u)")

hold off

%% plot

% uprime
uprime = avg_u_1_min - ubar_vec_array;

% plot(avg_time_1_min,ubar)

ubar_vector = zeros(length(u),1);
% replace each element by ubar value
ubar_vector(:,1) = ubar;
uprime = u - ubar_vector;

% bar u'u'
% vector u'u'
uprime_X_uprime = uprime.*uprime;
% mean u'u' bar over all time
bar_uprime_uprime = mean(uprime.*uprime);

% bar u'u' over 500 points
bar_uprime_uprime_500 = mean(uprime_X_uprime(1:500,1));

%% plot u'u' bar over running time interval 15 min

% % u' u' bar over first 15 min time interval
% n=1;
% for m=0:15*60*10:3
%     uprime_X_uprime_1min_avg(n,1) = mean(uprime_X_uprime(1,(1+m:m+15*60*10)))
%     n=n+1;
% end

% u'u' bar over first 15 min time interval
n=1;
for m=0:15*60*10:252000-(15*60*10)
    uprime_X_uprime_1min_avg(n,1) = mean(uprime_X_uprime((1+m:m+15*60*10),1));
    n=n+1;
end
figure(2)
% plot u'u' bar over running time interval 15 min
plot(avg_time_1_min,uprime_X_uprime_1min_avg)
title("u'u' bar Vs time running mean 15 min")
xlabel('Time average 15 min')
ylabel("bar(u'u')")

%% Scatter plot u v 500 points

figure(3)
u_scatter((1:500),1) = uprime([1:500],1);
v_scatter((1:500),1) = uprime([501:1000],1);
scatter(u_scatter,v_scatter)
title("u Vs v scatter plot")


%% tau vary
% tau: time shift spacing
del_t=[10:10:2000];
tau=(10:10:2000)';
% bar over number u points
bar_over = 252000-2000;      % change bar_over to get autocorrelation graphs for different time intervals

%bar_over = [500 1000 10000 100000 252000-500];
bar_uprime500_100000_x_uprime_del_t_500_100000=zeros(200,1);

    for j=1:length(del_t)
    
        % bar_over 500 1000... 100000 points of u'
        uprime500_100000 = uprime([1:bar_over(1)],1);
        uprime_del_t_500_100000 = uprime([1+del_t(j):bar_over(1)+del_t(j)],1);
    
        % bar  u'u'(del t)
        uprime500_100000_x_uprime_del_t_500_100000 = uprime500_100000.*uprime_del_t_500_100000;
        bar_uprime500_100000_x_uprime_del_t_500_100000(j,1) = mean(uprime500_100000_x_uprime_del_t_500_100000);
        
        % normalize - divide u'u' square
        uprime_200 = uprime(1:length(bar_uprime500_100000_x_uprime_del_t_500_100000),1);
        uprime_uprime_square = uprime_200.*uprime_200;
        normalize_bar_uprime500_100000_x_uprime_del_t_500_100000 = bar_uprime500_100000_x_uprime_del_t_500_100000./(uprime_uprime_square);
       
    end

    hold on
    %plot(uprime500_100000_x_uprime_del_t_500_100000)
    figure(4)
    plot(tau,bar_uprime500_100000_x_uprime_del_t_500_100000)
    title('Autocorrelation of R(tau) Vs Tau')
    xlabel('Tau')
    ylabel('R(Tau)')
    %plot(uprime500_100000)
    hold off

%% pdf curve
dx=.15;
bins=0:dx:max(u);
freq=zeros(1,length(bins));

% -1 so element not exceed than 1
for j=1:length(bins)-1
    count=0;
    for i=1:length(u)
        if u(i,1)>bins(j) && u(i,1)<bins(j+1)
            count=count+1;
            freq(j)=count;
        end  
    end
end

% divide by bins to cut oscillationsin graph - web
prob = freq./(252000*dx);
figure(5)
plot(bins,prob)
title('Proability density function (pdf) over time')
xlabel('Number of bins')
ylabel('Probability')

disp(freq);

%integral
int= trapz(bins,prob);

%mean

mean_int = trapz(bins,prob.*bins);
disp(mean_int);
%plot(bins)
