
%% Plot of distance-angle function to describe a boundary %%

K = imread('K.bmp');
rtheta(K);

function [r, theta] = rtheta(Iin)

    [y, x] = find(Iin);
    centroid = [mean(x), mean(y)];

    dx = x - centroid(1);
    dy = y - centroid(2);
    [theta, r] = cart2pol(dx, dy); % cartesian to polar

    theta = rad2deg(theta); % Convert theta to degrees
    theta = theta - min(theta); % start from 0

    % Plot the graph of r against theta
    figure;
    plot(theta, r, '.');
    title('Graph of r against \theta');
    xlabel('\theta (degrees)');
    ylabel('r');
    
    xlim([0 360])
    
end

