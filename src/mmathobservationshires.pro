;+
; :Description:
;   Saves visualisations of the hires NSO/Kitt Peak synoptic maps for
;   CR1913-CR1999 (1800x900)
;-
pro mmathObservationsHires
    ;+
    ; I/O Parameters
    ;-

    ;+ Colour table to use
    loadct, 71

    ;+ Directory containing FITS (to read)
    fit_dir = "./data/mag.hires/"

    ;+ Filepath of the Carrington Rotation metadata (to read)
    carr_file = "./data/dates/carrington.rotations"

    ;+ Root export directory
    exportdir = "./examples/magnetic-observations/"

    ;+ Directory to read/write script variables to
    variables_path = exportdir + "processed-field-hires.sav"


    ;+ Directory to save map visualisations to
    rotation_dir = exportdir + "maps-hires/"

    ;+
    ; End of parameters
    ;-

    ;+ Create output directories if they doesn't exist
    file_mkdir, exportdir
    file_mkdir, rotation_dir

    ;+ Extract rotation number from filenames (assumes they're correct)
    filenames = file_basename(file_search(fit_dir + '*.fits'))

    ;+ Recover Carrington rotation number from fit filenames
    rot = stregex(filenames, '[0-9]+', /extract)
    n_rot = n_elements(rot) ;number of carrington rotations

    ;+ Read other Carrington rotation metadata
    carrdata, rot, carr_file, date=date, yy=yy

    ;+ Prepare array to store each synoptic map
    data = dblarr(1800, 900, n_rot)

    ;+ Check if the fits have been saved to file already to avoid re-reading them
    filetest = file_test(variables_path)

    ;+ Check if fits have been processed already
    if (filetest eq 1) then begin
        ;+ Restore saved fit data
        restore, variables_path
    endif else begin
        ;+ Read in the fits and save image of each one
        temp = dblarr(1800, 900)
        for i=0, n_rot-1 do begin
            fxread, fit_dir+filenames[i], temp
            data[*, *, i] = temp
            savmap, temp, rotation_dir + string(rot[i]) + ".eps", 1.0, -10, 10
        endfor

        ;+ Save processed magnetic data for future reading
        save, data, filename=variables_path
    endelse
end
