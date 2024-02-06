FUNCTION RemoveRows, array, rows

  ; array -- A 2D array from which rows will be removed.
  ; rows -- A vector of row indices to remove from array.

  ; Need both positional parameters.
  IF N_Params() NE 2 THEN BEGIN
    Print, "Usage: 'RemoveRows, array, rowsToRemove'"
    RETURN, -1
  ENDIF

  ; The array must be 2D.
  ndims = Size(array, /N_DIMENSIONS)
  IF ndims NE 2 THEN BEGIN
    void = Dialog_Message('Array must be 2D.')
    Print, "Usage: 'RemoveRows, array, rowsToRemove'"
    RETURN, -1
  ENDIF

  ; The rows must be a vector.
  IF Size(rows, /N_DIMENSIONS) EQ 0 THEN rows = [rows]

  ; Find the dimensions of the array.
  dims = Size(array, /DIMENSIONS)

  ; Return the shortened array.
  RETURN, array[*, Where(~Histogram(rows, MIN=0, MAX=dims[1]-1), /NULL)]

