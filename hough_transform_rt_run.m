
%% Hough transform for line detection in the rho-theta transform space (normal representation) %%

letter = imread('letter.bmp');
letter_edge = edge(letter, 'canny', 0.5);

tic; % start timer

hough_transform_rt(letter_edge, 1, pi/250, 20) % rho-theta | output is saved as letter.bmp

elapsed_time = toc; % end timer
disp(['Time taken: ' num2str(elapsed_time) ' seconds']); % print time taken

% modified rho-theta with the line segments aligned with the object boudnary
% ht_plot_aligned(letter, letter_edge, 1, pi/250, 20) 

function hough_transform_rt(edge_image, rho_res, theta_res, n_lines)

    [height, width] = size(edge_image);
    max_rho = round(sqrt(height^2 + width^2));
    min_rho = -max_rho;
    rho_range = min_rho:rho_res:max_rho; % Range of rho values (rows)
    theta_range = -pi/2:theta_res:pi/2; % Range of theta values (columns)
    accumulator = zeros(length(rho_range), length(theta_range)); % Initialize the accumulator array with zeros

    [row, col] = find(edge_image); % Row and col index for all non-zero edge points (xy space)
    
    for i = 1:length(row) % Loop through these edge points
        y = row(i);
        x = col(i);
        
        for j = 1:length(theta_range) % for each theta
            rho = x * cos(theta_range(j)) + y * sin(theta_range(j)); % Calculate rho for the current theta
            
            % round to the nearest value
            for k = 1:length(rho_range) 
                closest_rho = abs(rho - rho_range(k)); % Find the closest matching rho value in the array
                if closest_rho <= 0.5 * rho_res
                    accumulator(k, j) = accumulator(k, j) + 1;
                end
            end
        end
    end

    regional_maximas = imregionalmax(accumulator);
    centers = findRegionCenters(regional_maximas);
    [maxima_rows, maxima_cols] = find(centers);
    maxima_triplets = [maxima_rows, maxima_cols, accumulator(sub2ind(size(accumulator), maxima_rows, maxima_cols))];
    sorted_maximas = sortrows(maxima_triplets, -3);
    top_maximas = sorted_maximas(1:min(n_lines, size(sorted_maximas, 1)), :);
    
    figure;
    imshow(edge_image);
    hold on;

    for k = 1 : size(top_maximas, 1)
        row = top_maximas(k, 1);
        col = top_maximas(k, 2);
        rho_value = rho_range(row);
        theta_value = theta_range(col);
        
        % Plot
        x_line = 1:width;
        y_line = (rho_value - x_line * cos(theta_value)) / sin(theta_value);
        plot(x_line, y_line, 'r', 'LineWidth', 2);

    end
    hold off;

end

function centers = findRegionCenters(regional_maximas)
    % Find the centers of all regions in a regional maximas array

    % Label the connected components (regions) in the regional maximas
    labeled_regions = bwlabel(regional_maximas);

    % Initialize an array of the same size as the input
    centers = zeros(size(regional_maximas));

    % Iterate through each labeled region
    for label = 1:max(labeled_regions(:))
        % Find the coordinates (row and column) of the region's pixels
        [row, col] = find(labeled_regions == label);

        % Calculate the centroid (center) of the region
        center_row = round(mean(row));
        center_col = round(mean(col));

        % Set the center location to 1 in the result array
        centers(center_row, center_col) = 1;
    end
end

function ht_plot_aligned(original_image, edge_image, rho_res, theta_res, n_lines)

    [height, width] = size(edge_image);
    max_rho = round(sqrt(height^2 + width^2));
    min_rho = -max_rho;
    rho_range = min_rho:rho_res:max_rho; % Range of rho values (rows)
    theta_range = -pi/2:theta_res:pi/2; % Range of theta values (columns)
    accumulator = zeros(length(rho_range), length(theta_range)); % Initialize the accumulator array with zeros

    [row, col] = find(edge_image); % Row and col index for all non-zero edge points (xy space)
    
    for i = 1:length(row) % Loop through these edge points
        y = row(i);
        x = col(i);
        
        for j = 1:length(theta_range) % for each theta
            rho = x * cos(theta_range(j)) + y * sin(theta_range(j)); % Calculate rho for the current theta
            
            % round to the nearest value
            for k = 1:length(rho_range) 
                closest_rho = abs(rho - rho_range(k)); % Find the closest matching rho value in the array
                if closest_rho <= 0.5 * rho_res
                    accumulator(k, j) = accumulator(k, j) + 1;
                end
            end
        end
    end

    regional_maximas = imregionalmax(accumulator);
    centers = findRegionCenters(regional_maximas);
    [maxima_rows, maxima_cols] = find(centers);
    maxima_triplets = [maxima_rows, maxima_cols, accumulator(sub2ind(size(accumulator), maxima_rows, maxima_cols))];
    sorted_maximas = sortrows(maxima_triplets, -3);
    top_maximas = sorted_maximas(1:min(n_lines, size(sorted_maximas, 1)), :);
    
    figure;
    imshow(original_image);
    hold on;

    for k = 1 : size(top_maximas, 1)
        tmp_row = top_maximas(k, 1);
        tmp_col = top_maximas(k, 2);
        rho_value = rho_range(tmp_row);
        theta_value = theta_range(tmp_col);

        % Plot
        x_line = 1:width;
        y_line = (rho_value - x_line * cos(theta_value)) / sin(theta_value);
        
        to_keep = true(size(x_line));
        % Compare y_line to the row indexes of edge_image
        for i = 1:length(x_line)
            idx = find(col == x_line(i));  % if x is in the (col value) vector of the edge point vector, take its index
            if ~isempty(idx) 
                corresponding_row = row(idx); % using the index, get the corresponding (row value)
                if ~any(abs(y_line(i) - corresponding_row) <= 2) % If this row value is not close to the corresponding value of y_line
                    to_keep(i) = false; % don't keep it when plotting
                end
            else
                to_keep(i) = false; % if x is not in the col value vector, don't keep that part of the time
            end
        end
        
        % All 1s before and after a long string of of 1s should become 0 
        % so that the line does not extend longer than it needs to
        to_keep = conv(double(to_keep), ones(1, 10), 'same') >= 5;

        plot(x_line(to_keep), y_line(to_keep), 'r', 'LineWidth', 2);
        % plot(x_line, y_line, 'r', 'LineWidth', 2);

    end
    hold off;

end
    
