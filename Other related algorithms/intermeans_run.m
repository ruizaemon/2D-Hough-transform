
%% Intermeans algorithm for automatic thresholding %%

I = imread('I.bmp');
intermeans(I);

function [T, Iout] = intermeans(Iin)

    T = mean(Iin(:));
    maxIterations = 100;

    for iter = 1:maxIterations
        % Partition the image into two groups, R1 and R2, using the current threshold T
        R1 = Iin(Iin <= T);
        R2 = Iin(Iin > T);
        
        mu_1 = mean(R1);
        mu_2 = mean(R2);

        newT = 0.5 * (mu_1 + mu_2);

        if abs(newT - T) < 0.01
            break;
        end

        T = newT;
    end

    Iout = Iin > T;
    imwrite(Iout, 'I1.bmp');

    % Plot histogram
    figure;
    subplot(2,2,1);
    imshow(Iin);
    title('Original image')

    subplot(2,2,2);
    imhist(Iin);
    h = imhist(Iin);
    ylim([0, max(h)*1.01]); % Set y-axis limits to cover the full range of frequencies
    hold on;
    plot([T, T], [0, max(h)], 'r', 'LineWidth', 2); % Add a red line at the final threshold
    hold off;
    title(['Histogram with intermeans T = ', num2str(T)]);
    
    subplot(2,2,3);
    imshow(Iout);
    title(['Thresholded image at T=', num2str(T)]);

end