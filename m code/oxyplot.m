function [h,ax]=oxyplot(time,d18o,varargin)
%OXYPLOT Plot an oxygen isotope record quickly
%
%SYNTAX:   
%       oxyplot(time,d18o)
%       where time is your time/age vector of increasing values
%       and d18O is your oxygen isotope data
%       Optional: add axis to plot onto
%
%       [h,ax]=oxyplot(time,d18o)
%       gives you h, the line handle and ax the axis handle

%preamble
x=time;
y=d18o;
if nargin==2
[h]=plot(x,y,'color','k');
end
if nargin==3  
    ax=varargin{1};
[h]=plot(ax,x,y,'color','k');
end

ax=gca;
set(ax,'ydir','reverse',...
    'xlim',[min(x) max(x)]);
set(gcf,'color','w');
box off
xlabel('Time (ka)');
ylabel(['\delta ^1^8O (',char(8240),')']); 

end

  