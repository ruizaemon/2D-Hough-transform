
%% Chain code to represent a boundary %%

K = imread('K.bmp');
ChainCode(K);

function [my_code] = ChainCode(Iin)

    [y, x] = find(Iin);
    
    % Sort the points based on polar angle in clockwise order
    centroid = [mean(x), mean(y)];
    angles = atan2(y - centroid(2), x - centroid(1));
    [~, order] = sort(angles);
    x = x(order);
    y = y(order);

    num_points = length(x);
    num_segments = 40; % manually defined
    segment_length = floor(num_points / num_segments);
    segment_endpoints_x = zeros(num_segments, 2);
    segment_endpoints_y = zeros(num_segments, 2);

    for i = 1:num_segments
        start_idx = (i - 1) * segment_length + 1;
        end_idx = i * segment_length;
        
        % Check for the last segment
        if i == num_segments
            end_idx = num_points;
        end

        % Assign segment endpoints
        segment_endpoints_x(i, 1) = x(start_idx);
        segment_endpoints_y(i, 1) = y(start_idx);
        segment_endpoints_x(i, 2) = x(end_idx);
        segment_endpoints_y(i, 2) = y(end_idx);
    end

    % Initialize variables to store the chain code
    my_code = '';

    figure;
    hold on;

    % Calculate the angle and assign the direction for each line segment
    for i = 1:num_segments
        % Calculate the angle of the line relative to east
        angle = atan2(segment_endpoints_y(i, 2) - segment_endpoints_y(i, 1), ...
                      segment_endpoints_x(i, 2) - segment_endpoints_x(i, 1));

        % Convert the angle to degrees
        angle_degrees = rad2deg(angle);

        % Map the angle to the corresponding direction
        direction = round(mod(angle_degrees + 360, 360) / 45);
        if direction == 8
            direction = 0;
        end

        % Append the direction to the chain code
        my_code = [my_code, num2str(direction)];

        % Plot the line segment
        plot([segment_endpoints_x(i, 1), segment_endpoints_x(i, 2)], ...
             [segment_endpoints_y(i, 1), segment_endpoints_y(i, 2)], '-o');

    end

    % Plotting
    axis equal;
    xlabel('X');
    ylabel('Y');
    title('Line Segments');
    hold off;

    % Starting point normalization
    min_code = my_code;
    for i = 1:length(my_code)
        shifted_code = circshift(my_code, [0, -i]); % Circular shift the chain code
        shifted_number = str2double(shifted_code);
        min_number = str2double(min_code);
        if shifted_number < min_number
            min_code = shifted_code;
        end
    end
    my_code = min_code % display code

end