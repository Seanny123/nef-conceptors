function par = mcplotframe_testDiscard(d, n, p, proj)
% Plots frames of motion capture data.
%
% syntax
% par = mcplotframe(d, n);
% par = mcplotframe(d, n, p);
% par = mcplotframe(d, n, p, proj);
%
% input parameters
% d: MoCap data structure
% n: vector containing the numbers of the frames to be plotted
% p: animpar structure
% proj: projection used: 0 = orthographic (default), 1 = perspective
%
% output
% par: animpar structure used for plotting the frames (if color strings were used, they will converted to RGB triplets)
%
% examples
% par = mcplotframe(d, 1);
% mcplotframe(d, 500:10:600, par, 1);
%
% comments
% If the animpar structure is not given, the function calls
% mcinitanimpar and sets the .limits field of the animpar structure
% automatically so that all the markers fit into all frames.
% If the par.pers field (perspective projection) is not given, it is created internally for backwards
% compatibility. For explanation of the par.pers field, see help mcinitanimpar
%
% see also
% mcanimate, mcinitanimpar
%
% © Part of the Motion Capture Toolbox, Copyright ©2008,
% University of Jyvaskyla, Finland

if nargin<4
    proj=0;
end

if nargin<3
    p = mcinitanimpar;
    if nargin==1 || ~isnumeric(n) %output fix [BB20111031]
        disp([10, 'Error: please set frame number(s) you want to plot.', 10]);
        return
    end
end

for k=1:length(n)
    n1=n(k);
    if n1>d.nFrames
        w1=sprintf('Error: indicated frame(s) (%d) exceeds number of frames in data (%d).', max(n), d.nFrames);
        disp([10, w1, 10]);
        return
    end
end




%for compatibility
if ~isfield(p, 'pers')
    p.pers.c=[0 -3000 0];
    p.pers.th=[0 0 0];
    p.pers.e=[0 -2000 0];
end

par=p;

% color management (compatibility): convert old string color definition into num array [BB20111031]
if ischar(p.colors)
    colors=NaN(5,3);
    for k=1:5
        colors(k,:)=lookup_l(p.colors(k));
    end
end

bgcol=colors(1,:); % set background color

if isfield(p,'markercolors') && ~isempty(p.markercolors) %field and colors are set
    if ischar(p.markercolors) % but in string format
        mcol=NaN(d.nMarkers,3);
        for k=1:size(p.markercolors,2)
            mcol(k,:)=lookup_l(p.markercolors(k));%convert to num array
        end
    else 
        mcol=p.markercolors; %field and colors are set in num format already
    end
    if d.nMarkers > length(p.markercolors)
        for k=length(p.markercolors)+1:d.nMarkers
            mcol(k,:)=colors(2,:);
        end
    end
else %no field or/and empty
    mcol=repmat(colors(2,:),d.nMarkers,1);
end

if isfield(p,'conncolors') && ~isempty(p.conncolors) %field and colors are set
    if ischar(p.conncolors) % but in string format
        ccol=NaN(size(p.conn,1),3);
        for k=1:size(p.conncolors,2)
            ccol(k,:)=lookup_l(p.conncolors(k));%convert to num array
        end
    else 
        ccol=p.conncolors; %field and colors are set in num format already
    end
    if size(p.conn,1) > length(p.conncolors)
        for k=length(p.conncolors)+1:size(p.conn,1)
            ccol(k,:)=colors(3,:);
        end
    end
else %no field or/and empty
    ccol=repmat(colors(3,:),size(p.conn,1),1);
end

%BBFIX 20120404: mcmerge problems
if p.trl~=0
    tcol=repmat(colors(4,:),d.nMarkers,1);
    if isfield(p,'tracecolors') && ~isempty(p.tracecolors) %(field and) tracecolors are set
        if ischar(p.tracecolors) % but in string format
            for k=1:size(p.tracecolors,2)
                tcol(k,:)=lookup_l(p.tracecolors(k));%convert to num array
            end
        else
            tcol(1:size(p.tracecolors,1),:)=p.tracecolors; %tracecolors in num format already
        end
    end
    if isempty(p.trm)
        p.trm=1:d.nMarkers; %plot all traces if trm is empty
    else
        if length(p.trm)<d.nMarkers
            %sort p.tracecolors/tcol1 according to p.trm
            %increase p.trm to have same length as d.nMarkers and fill it up with NaNs
            tmp=repmat(colors(4,:),d.nMarkers,1);
            tmp1=nan(d.nMarkers,1)';
            for k=1:length(p.trm)
                tmp(p.trm(k),:)=tcol(k,:);
                tmp1(p.trm(k))=p.trm(k);
            end
            tcol=tmp;
            p.trm=tmp1;
        end
    end
else
    tcol=[];
    p.trm=[];
end
p.tracecolors=tcol;

if ~isempty(p.trm)%created problem when merging and translating...
    p.trm=p.trm(:);
end

%warning if trace length is empty, but marker trace vector is set 
if (p.trl==0 && ~isempty(p.trm)) || (p.trl==0 && ~isempty(p.tracecolors))
    disp([10, 'Warning: Please set trace length (trl) in your animation parameters in order to plot traces.', 10])
end

%if trace length is set, but vector indicating the markers to be traced is empty, all markers will be trace
if p.trl~=0 && isempty(p.trm)
    disp([10, 'Note: All markers traced.', 10])
    p.trm=1:d.nMarkers;
end


%BBFIX20120404: mcmerge problems
if p.showmnum==1
    ncol=repmat(colors(5,:),d.nMarkers,1);
    if isfield(p,'numbercolors') && ~isempty(p.numbercolors) %(field and) numbercolors are set
        if ischar(p.numbercolors) % but in string format
            for k=1:size(p.numbercolors,2)
                ncol(k,:)=lookup_l(p.numbercolors(k));%convert to num array
            end
        else
            ncol(1:size(p.numbercolors,1),:)=p.numbercolors; %numbercolors in num format already
        end
    end
    if isempty(p.numbers)
        p.numbers=1:d.nMarkers; %plot all markers if numbers is empty
    else
        if length(p.numbers)<d.nMarkers
            %sort p1.numbercolors/ncol1 according to p1.numbers
            %increase p1.numbers to have same length as d1.nMarkers and fill it up with NaNs
            tmp=repmat(colors(5,:),d.nMarkers,1);
            tmp1=nan(d.nMarkers,1)';
            for k=1:length(p.numbers)
                tmp(p.numbers(k),:)=ncol(k,:);
                tmp1(p.numbers(k))=p.numbers(k);
            end
            ncol=tmp;
            p.numbers=tmp1;
        end
    end
else
    ncol=[];
    p.numbers=[]; %fill up p1.numbers with nan
end


par.colors=colors;
par.markercolors=mcol;
par.conncolors=ccol;
par.tracecolors=tcol;
par.numbercolors=ncol;

par.numbers=p.numbers;


%fill up trace widths (twidth) to have same length as traced markers (trm) / nMarkers
if p.trl~=0
    if length(p.twidth)<d.nMarkers
        twidth=nan(d.nMarkers,1)';
        i=1;
        for k=1:length(p.trm)
            if isnan(p.trm(k))
            else
                twidth(k)=p.twidth(i);
                if i<length(p.twidth)
                    i=i+1;
                end
            end
        end
        p.twidth=twidth;
    end
end
    
%individual widths for connectors and traces
%fill up cwidth
if length(p.cwidth)<size(p.conn,1)
    cl=length(p.cwidth);
    cwidth=nan(size(p.conn,1),1);
    cwidth(1:length(p.cwidth))=p.cwidth;
    for k=cl+1:length(cwidth)
        if isnan(cwidth(k))
            cwidth(k)=cwidth(cl);%fill up with last given value
        end
    end
elseif length(p.cwidth)>size(p.conn,1)%if cwidth is longer than conn
    cwidth=p.cwidth(1:size(p.conn,1));
else
    cwidth=p.cwidth;
end
p.cwidth=cwidth;

% if isfield(p,'cwidth') && length(p.cwidth)>1
%     p.cwidth=p.cwidth;
% else
%     p.cwidth=repmat(p.cwidth(1),1,length(p.conn));
% end
if isfield(p,'twidth') && length(p.twidth)>1
    p.twidth=p.twidth;
else
    if ~isempty(p.trm)
        p.twidth=repmat(p.twidth(1),1,length(p.trm));
    end
end
 


az=p.az;
el=p.el;
    
d1 = mcrotate(d, az, [0 0 1]);
d2 = mcrotate(d1, el, [1 0 0]); %%%%%%%
    
if proj==0 % orthographic projection    
    x=d2.data(n,1:3:end);
    y=d2.data(n,2:3:end);
    z=d2.data(n,3:3:end);
else % perspective projection
    if ~isfield(p,'pers') % for backward compatibility, use default values
        p.pers.c=[0 -4000 0];
        p.pers.th=[0 0 0];
        p.pers.e=[0 -2000 0];
    end

    th=180*p.pers.th/pi;
    rot1=[1 0 0; 0 cos(th(1)) -sin(th(1)); 0 sin(th(1)) cos(th(1))];
    rot2=[cos(th(2)) 0 sin(th(2)); 0 1 0; -sin(th(2)) 0 cos(th(2))];
    rot3=[cos(th(3)) -sin(th(3)) 0; sin(th(3)) cos(th(3)) 0; 0 0 1];
    dd=zeros(size(d2.data(n,:)));
    %p.pers.e(2)=min(min(d.data(n,2:3:end)));
    for k=1:d.nMarkers
        dd(:,3*k+(-2:0))=(rot1*rot2*rot3*(d2.data(n,3*k+(-2:0))'-repmat(p.pers.c',1,length(n))))';
    end
    % make closest marker to be on the projection plan
    dd(:,2:3:end)=dd(:,2:3:end)-min(min(dd(:,2:3:end)))+p.pers.e(2)-p.pers.c(2); 
    x=-(dd(:,1:3:end)-repmat(p.pers.e(1),length(n),d.nMarkers)).*(p.pers.e(2)./dd(:,2:3:end));
    y=(p.pers.e(2)-p.pers.c(2))./dd(:,2:3:end); % used for marker size scaling
    z=-(dd(:,3:3:end)-repmat(p.pers.e(3),length(n),d.nMarkers)).*(p.pers.e(2)./dd(:,2:3:end));    
end


if isempty(p.limits)
    % find ranges of coordinates
    tmp=x(:); maxx=max(tmp(find(~isnan(tmp)))); minx=min(tmp(find(~isnan(tmp))));
    tmp=y(:); maxy=max(tmp(find(~isnan(tmp)))); miny=min(tmp(find(~isnan(tmp))));
    tmp=z(:); maxz=max(tmp(find(~isnan(tmp)))); minz=min(tmp(find(~isnan(tmp))));
    midx = (maxx+minx)/2;
    midy = (maxy+miny)/2;
    midz = (maxz+minz)/2;
    
    scrratio = p.scrsize(1)/p.scrsize(2);
    range = max((maxx-minx)/scrratio, maxz-minz)/2;
    zrange = (maxy-miny)/2;
    % axis limits for plot
    p.limits = [midx-scrratio*1.2*range midx+scrratio*1.2*range midz-1.2*range midz+1.2*range];
end
minxx = p.limits(1);
maxxx = p.limits(2);
minzz = p.limits(3);
maxzz = p.limits(4);

for k=1:size(x,1) % main loop
    if 0
        clf
    else
        figure;
    end
    set(gcf,'Position',[50 50 p.scrsize(1) p.scrsize(2)]) ; % DVD: w=720 h=420
    
    axes('position', [0 0 1 1]); hold on
    set(gcf, 'color', bgcol);
    view(0,90);
    colormap([ones(64,1) zeros(64,1) zeros(64,1)]);
    
    % plot marker-to-marker connections
    if ~isempty(p.conn)
        for m=1:size(p.conn,1)
            %if x(k,p.conn(m,1))*x(k,p.conn(m,2))~=0
            if isfinite(x(k,p.conn(m,1))*x(k,p.conn(m,2)))
                plot([x(k,p.conn(m,1)) x(k,p.conn(m,2))], [z(k,p.conn(m,1)) z(k,p.conn(m,2))], '-','Color',ccol(m,:),'LineWidth', p.cwidth(m));
            end
        end
    end
    
    % plot midpoint-to-midpoint connections
    if ~isempty(p.conn2)
        for m=1:size(p.conn2,1)
            %if x(k,p.conn2(m,1))*x(k,p.conn2(m,2))*x(k,p.conn2(m,3))*x(k,p.conn2(m,4))~=0
            if isfinite(x(k,p.conn2(m,1))*x(k,p.conn2(m,2))*x(k,p.conn2(m,3))*x(k,p.conn2(m,4)))
                tmpx1 = (x(k,p.conn2(m,1))+x(k,p.conn2(m,2)))/2;
                tmpx2 = (x(k,p.conn2(m,3))+x(k,p.conn2(m,4)))/2;
                tmpy1 = (z(k,p.conn2(m,1))+z(k,p.conn2(m,2)))/2;
                tmpy2 = (z(k,p.conn2(m,3))+z(k,p.conn2(m,4)))/2;
                plot([tmpx1 tmpx2], [tmpy1 tmpy2], '-','Color',ccol(m,:),'LineWidth', p.cwidth(m));
            end
        end
    end
    
    % plot traces if animation
    
    
    
    % plot markers
    for m=1:size(x,2)
        %if x(k,m)~=0 & ~isnan(x(k,m)) % if marker visible
        if isfinite(x(k,m)) % if marker visible
            if proj==0 % orthographic projection
                plot(x(k,m),z(k,m),'o','MarkerSize',p.msize(min(m,length(p.msize))),'MarkerEdgeColor',mcol(m,:),'MarkerFaceColor',mcol(m,:))
            else % perspective projection
%                plot(x(k,m),z(k,m),[mcol(m) 'o'],'MarkerSize',p.msize(min(m,length(p.msize))),'MarkerFaceColor',mcol(m))
                plot(x(k,m),z(k,m),'o','MarkerSize',round(y(k,m)*p.msize(min(m,length(p.msize)))),'MarkerEdgeColor',mcol(m,:),'MarkerFaceColor',mcol(m,:))
            end
            if p.showmnum
                if isempty(p.numbers)
                    h=text(x(k,m),z(k,m)+50,num2str(m));
                    set(h,'FontSize',16);
                    set(h,'Color',ncol(m,:))
                else
                     %if ismember(m, p.numbers) %FIX BB20120326: redefinining numbers plotting (mcmerge) 
                     if m<=length(p.numbers)
                         if isnan(p.numbers(m)) %NaN numbers not plotted
                         else h=text(x(k,m),z(k,m)+50,num2str(p.numbers(m)));
                             set(h,'FontSize',16);
                             set(h,'Color',ncol(m,:))
                         end
                    end
                end
            end
            axis off
        end
        axis([minxx maxxx minzz maxzz]); axis off
    end
    if p.showfnum
        h=text(minxx+0.95*(maxxx-minxx), minzz+0.97*(maxzz-minzz), num2str(k),...
            'HorizontalAlignment','Right','FontSize',12,'FontWeight','bold');
        set(h,'Color',colors(5,:))
    end
    
%     %plot some copywrite text or anything else - BB20121102
%     text(maxxx-250, 0, '©', 'FontSize', 40, 'color', [.7 .7 .7]);
%     text(maxxx-300, 0, {'Birgitta Burger', 'Jyväskylä Univ.', 'Finland'});
%     text(minxx+40, minzz+0.97*(maxzz-minzz), 'High Sub-Band 2 Flux', 'FontSize', 16, 'FontWeight', 'bold');
%     text(minxx+70, minzz+0.97*(maxzz-minzz)-75, {'high speed of head'}, 'FontSize', 12, 'FontWeight', 'bold');
    
    drawnow
    hold off
    
      
end



return;


function colorar=lookup_l(colorstr)
if strcmp(colorstr, 'k')
    colorar=[0 0 0];
elseif strcmp(colorstr, 'w')
    colorar=[1 1 1];
elseif strcmp(colorstr, 'r')
    colorar=[1 0 0];
elseif strcmp(colorstr, 'g')
    colorar=[0 1 0];
elseif strcmp(colorstr, 'b')
    colorar=[0 0 1];
elseif strcmp(colorstr, 'y')
    colorar=[1 1 0];
elseif strcmp(colorstr, 'm')
    colorar=[1 0 1];
elseif strcmp(colorstr, 'c')
    colorar=[0 1 1];
end

return;



