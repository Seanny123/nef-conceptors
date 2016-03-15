function mcplotphaseplane(d1, d2, p1, p2)
% Plots motion capture data on a phase plane.
%
% syntax
% mcplotphaseplane(d1, d2, marker, dim) % for MoCap data structure
% mcplotphaseplane(n1, n2, marker) % for norm data structure
% mcplotphaseplane(s1, s2, segm, var) % for segm data structure
%
% input parameters
% d1, d2, n1, n2, s1, s2: MoCap data structure, norm data structure, or segm data structure
% marker: vector containing marker numbers to be plotted (for MoCap and norm data structure)
% dim: vector containing dimensions to be plotted (for MoCap data structure)
% segm: body segment number (for segm data structure)
% var: variable to be plotted for segment segm (for segm data structure)
%
% output
% Figure.
%
% examples
% mcplotphaseplane(d1, d2, 1:3, 3) % for MoCap data structure
% mcplotphaseplane(n1, n2, 5) % for norm data structure
% mcplotphaseplane(s1, s2, [2 6 25], 'angle') % for segm data structure
% mcplotphaseplane(s1, s2, 5:10, 'eucl') % for segm data structure
% mcplotphaseplane(s1, s2, [12 14], 'quat') % for segm data structure
%
% © Part of the Motion Capture Toolbox, Copyright ©2008,
% University of Jyvaskyla, Finland

if isfield(d1,'type') && isfield(d2, 'type')
    if strcmp(d1.type, d2.type) && d1.freq==d2.freq && d1.nFrames==d2.nFrames
        t = (1:d1.nFrames)'-1;
        if strcmp(d1.type, 'MoCap data')
            for k=1:length(p1)
                for m=1:length(p2)
                    subplot(length(p1), length(p2), length(p2)*(k-1)+m)
                    plot(d1.data(:,3*p1(k)-3+p2(m)), d2.data(:,3*p1(k)-3+p2(m)))
                    axis([-Inf Inf -Inf Inf])
                    title(['Marker ' num2str(p1(k)) ', dim.' num2str(p2(m))])
                end
            end
        elseif strcmp(d1.type, 'norm data')
            plot(d1.data(:,p1), d2.data(:,p1));
            for k=1:length(p1)
                subplot(length(p1), 1, k)
                plot(d1.data(:,p1(k)), d2.data(:,p1(k)))
                axis([-Inf Inf -Inf Inf])
                title(['Marker ' num2str(p1(k))])
            end
        elseif strcmp(d1.type, 'segm data')
            tmp=[];
            for k=1:length(p1)
                tmp1 = getfield(d1.segm(p1(k)),p2);
                tmp2 = getfield(d2.segm(p1(k)),p2);
                if ~isempty(tmp1)
                    if strcmp(p2, 'angle')
                        subplot(length(p1),1,k)
                        plot(tmp1, tmp2)
                        axis([-Inf Inf -Inf Inf])
                        title(['Segm. ' num2str(p1(k))])
                    elseif strcmp(p2, 'eucl')
                        for m=1:3
                            subplot(length(p1), 3, 3*(k-1)+m)
                            plot(tmp1(:,m), tmp2(:,m)), title(['Segm. ' num2str(p1(k)) ' Dim. ' num2str(m)]);
                            axis([-Inf Inf -Inf Inf])
                        end
                    elseif strcmp(p2, 'quat')
                        for m=1:4
                            subplot(length(p1), 4, 4*(k-1)+m)
                            plot(tmp1(:,m), tmp2(:,m)), title(['Segm. ' num2str(p1(k)) ' Comp. ' num2str(m)]);
                            axis([-Inf Inf -Inf Inf])
                        end
                    end
                else
                    disp('No data to be plotted.')
                end
            end
        end
    end
else % direct reference to data
    if ~exist('p1')
        p1=1;
    end
    plot(d1(:,p1),d2(:,p1));
end

