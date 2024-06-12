
%% Hough transform for line detection in the ab (slope-intercept) transform space %%

letter = imread('letter.bmp');
letter_edge = edge(letter, 'canny', 0.5); % compute edge map

tic; % start timer

hough_transform_ab(letter_edge, 0.2, 1, 20) % slope-intercept

elapsed_time = toc; % end timer
disp(['Time taken: ' num2str(elapsed_time) ' seconds']); % print time taken

function hough_transform_ab(edge_image, a_res, b_res, n_lines)
    % Size of the input image
    [height, width] = size(edge_image);
    
    % Range of a values
    max_a = round(sqrt(height^2 + width^2));
    min_a = -max_a/2;
    a_range = min_a:a_res:max_a;
    
    % Range of b values
    max_b = max_a*4;
    min_b = -max_a;
    b_range = min_b:b_res:max_b;
    
    % Initialize the accumulator array with zeros
    accumulator = zeros(length(a_range), length(b_range));
    
    % Row and col index for all non-zero edge points (xy space)
    [row, col] = find(edge_image);
    
    % Loop through these edge points
    for i = 1:length(row)
        y = row(i);
        x = col(i);
        b = -1*x*a_range + y;  % b = -x*a+y
        
        % Loop through each b value and round to the nearest value
        for j = 1:length(b)
            
            % Find the closest matching b value in the array
            [closest_b, b_index] = min(abs(b(j) - b_range));
            
            if closest_b <= 0.5 * b_res
                % Increment the corresponding cell in the accumulator
                accumulator(j, b_index) = accumulator(j, b_index) + 1;
            end
        end
    end

    % Find regional maximas in the accumulator
    regional_maximas = imregionalmax(accumulator);

    % Find the centers of these regional maximas
    centers = findRegionCenters(regional_maximas);

    % Get the indices of regional maximas
    [maxima_rows, maxima_cols] = find(centers);
    
    % Create an array of (row, col, value) triplets
    maxima_triplets = [maxima_rows, maxima_cols, accumulator(sub2ind(size(accumulator), maxima_rows, maxima_cols))];

    % Sort the regional maximas by their values in descending order
    sorted_maximas = sortrows(maxima_triplets, -3);

    % Take the top n_lines maximas
    top_maximas = sorted_maximas(1:min(n_lines, size(sorted_maximas, 1)), :);

    figure;
    imshow(edge_image); % Display the edge image
    hold on;

    for k = 1 : size(top_maximas, 1)
        row = top_maximas(k, 1);
        col = top_maximas(k, 2);

        % Get the corresponding 'a' and 'b' values from the indices
        a_value = a_range(row);
        b_value = b_range(col);

        % Plot the detected line on the image
        x_line = 1:width;
        y_line = a_value * x_line + b_value;
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