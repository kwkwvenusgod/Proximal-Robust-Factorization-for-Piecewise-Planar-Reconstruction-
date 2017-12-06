function [ translation_error, translation_error_average, rotation_angular_error,rotation_angular_error_average, rotation_epe, rotation_epe_average] = Evaluation_Motion( gt, est )
    %EVALUATION_MOTION Summary of this function goes here
    %   Detailed explanation goes here
    % Normalize
    gt_translation = gt(:,1:3);
    gt_translation_mag = sqrt(sum(gt_translation.*gt_translation,2));
    gt_translation_normalize = gt_translation./repmat(gt_translation_mag,1,3);

    est_translation = est(:,1:3);
    est_translation_mag = sqrt(sum(est_translation.*est_translation,2));
    est_translation_normalize = est_translation./repmat(est_translation_mag,1,3);

    translation_error = sum(gt_translation_normalize.*est_translation_normalize,2);
    translation_error = acos(translation_error);
    translation_error = 180 * translation_error/pi;
    translation_error_average = mean(translation_error);

    gt_rotation = gt(:,4:end);
    gt_rotation_mag = sqrt(sum(gt_rotation.*gt_rotation,2));
    gt_rotation_normalize = gt_rotation./repmat(gt_rotation_mag,1,3);

    est_rotation = est(:,4:end);
    est_rotation_mag = sqrt(sum(est_rotation.*est_rotation,2));
    est_rotation_normalize = est_rotation./repmat(est_rotation_mag,1,3);

    rotation_angular_error = sum(gt_rotation_normalize.*est_rotation_normalize,2);
    rotation_angular_error = acos(rotation_angular_error);
    rotation_angular_error = 180 *rotation_angular_error /pi;
    rotation_angular_error_average = mean(rotation_angular_error);
    
    
    rotation_epe = gt_rotation - est_rotation;
    rotation_epe = rotation_epe.*rotation_epe;
    rotation_epe = sum(rotation_epe,2);
    rotation_epe = sqrt(rotation_epe);
    rotation_epe_average = mean(rotation_epe);
end

