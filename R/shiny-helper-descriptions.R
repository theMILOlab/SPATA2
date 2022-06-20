
# createImageAnnotations()
create_image_annotations_descr <- list(

  caption = c("Display the image annotation tags as a caption. (If display mode: One by one.)"),
  encircle = c("Display the polygon with which the structure has been encircled. (If display mode: One by one.)"),
  expand = c(""),
  color_by = c("Use SPATA variables to color the surface of the image."),
  display_mode =
    c(
      "If 'Surface', the image annotations are projected on the whole histology image.",
      "",
      "If 'One by one', each image annotations is displayed in a separate window."
    ),
  drawing_mode =
    c(
      "Single: Mode that allows highly specific encircling. At any time you can stop the drawing by exiting
    drawing via keyboard-shortcut 'e' or double clicking on the plot to adjust zooming. To continue drawing
    enter drawing via keyboard-shortcut 'd' or double clicking on the plot. The end of the previous line
    will connect to the current position of your cursor and you can continue encircling the structure.
    Additionally, you can provide a specific image annotation ID.",
      "",
      "Multiple: Mode that allows to quickly encircle several similar structures, thus adding multiple image
    annotations at the same time. Exiting the drawing mode immediately highlights the structure and
    entering it again starts the encircling of a new structure. Image annotation IDs are created as
    a combination of 'img_ann' and the position the annotations have in the list of image annotations."
    ),
  img_ann_id = c("The ID that uniquely identifes the image annotation."),
  img_ann_ids_select = c("The image annotations that you want to include in the plot above."),
  img_ann_tags =
    c(
      "Tag image annotations with specific words to describe their characteristics based on
      which they can be filtered or sorted."
    ),
  img_ann_tags_select = c(
    "Use tags to filter the image annotations below. Depending on the tags chosen and
      the option to handle them (select option on the left next to the tag list) the
      image annotation IDs below are filtered. There are four options to handle the
      image annoation filtering:",
    "",
    "1. 'ignore': All image annoations are selectable irrespective of the selected tags.",
    "",
    "2. 'any': Image annotations are selectable if they contain at least one
      of the chosen tags.",
    "",
    "3. 'all': Image annotations are selectable if they contain all of the chosen tags.",
    "",
    "4. 'identical': Image annotations are selectable if they match the chosen tags exactly."
  ),
  linesize = c("The size of the line that encircles the annotated structure."),
  ncol = c("Number of cols in which the windows are displayed. Ignored if 0."),
  nrow = c("Number of rows in which the windows are displyed. Ignored if 0."),
  pick_action_single =
    c(
      "Highlight: Closes the drawn circle and highlights the area that it includes marking
    the exact structure that you have annotated. Before you cann add an image annotation
    it must be higlighted. Alternatively you can use the keyboard-shortcut 'h'.",
      "",
      "Reset: Removes any drawing that is currently displayed on the interactive plot.
    Alternatively you can use the keyboard-shortcut 'r'."
    ),
  pick_action_multiple =
    c(
      "Reset all: Removes any drawing that is currently displayed on the interactive plot
    including already highlighted structures. Alternatively you can use the keyboard-shortcut 'a'.",
      "",
      "Reset last: Removes the most recently highlighted annotation. Alternatively you can
    use the keyboard-shortcut 'l'."
    ),
  pointsize = c("The size with which points are displayed if color the surface by a SPATA2 variable."),
  square = c("Force the image annoatation to be displayed in a square. (If display mode: One by one.)"),
  subtitle = c("Display the image annotation ID. (If display mode: One by one.)"),
  tab_panel_image_annotations = c(
    "Overview about the image annotations that are currently stored in the SPATA2 object."
  ),
  tab_panel_interaction =
    c(
      "This plot allows to interactively zoom in and out on your sample as well as to encircle the regions you want to name.",
      "",
      "Zooming: Brush the area on the plot you want to zoom in on. Then click on 'Zoom in'. You can zoom stepwise. To zoom one step
       back click on 'Zoom back'. To zoom out completely click on 'Zoom out'. Note that you can not zoom if
      you are drawing. If you want to stop drawing to zoom in on the image exit the drawing mode via shortcut 'e'
      zoom in and then start drawing again via doubleclicking or pressing 'd'. This is only possible if
      you are using drawing mode 'Single'.",
      "",
      "Encircling: By doubleckling or pressing 'd' you start drawing. Encircle the area/structure you want to annotate by simply moving
      the cursor. By double clicking again or pressing 'e' you stop drawing. Depending on the drawing mode you have chosen (Single
      or Multiple) the encircled area is highlighted immediately (Multiple) or you need to click on 'Highlight' or press 'h' (Single).",
      "",
      "Tagging: Provide additional information about the annotated structure in form of bullet points that can be used later on
      to group and/or separate them.",
      "",
      "Naming: This creates the image annotation ID to uniquely identify each annotation. If you chose the drawing mode 'Single' you can name the annotated structure individually. If you are using 'Multiple' the
      names are automatically generated as a combination of 'img_ann' and a number."
    ),
  tab_panel_orientation = c(" This plot mainly stays as is. Once you start zooming in on the interactive plot a rectangle is drawn to visualize where you currently are."),
  title = c("Display the number of the image annotation. (If display mode: One by one.)"),
  transparency = c("Change the transparency of the surface with which the annotated structure is highlighted."),
  transparency_point = c("The transparency of the points if you color the surface by a SPATA2 variable."),
  zooming_options = c(
    "Brush the area on the plot you want to zoom in on. Then click on 'Zoom in'. You can zoom stepwise. To zoom one step
    back click on 'Zoom back'. To zoom out completely click on 'Zoom out'. Note that you can not zoom if
    you are drawing. If you want to stop drawing to zoom in on the image exit the drawing mode via shortcut 'e'
    zoom in and then start drawing again via doubleclicking or pressing 'd'. This is only possible if
    you are using drawing mode 'Single'."
  )
)


# createSegmentation()

create_segmentation_descr <- list(

  color_by = c("Use SPATA variables to color the surface of the image."),
  linesize = create_image_annotations_descr$linesize,
  pick_action_interaction =     c(
    "Highlight: Closes the drawn circle and highlights the area that it includes marking
    the exact structure that you have annotated. Before you cann add an image annotation
    it must be higlighted.",
    "",
    "Reset: Removes any drawing that is currently displayed on the interactive plot."
  ),

  pick_action_overview = c(
  "The segment chosen on the left under 'Choose a group/segment'. Can either be renamed
   or discarded. Clicking on either of the two buttons opens a model in which to specify
   the action."),
  plot_interaction = c(
    "This plot allows to interactively zoom in and out on your sample as well as to encircle the regions you want to name.",
    "",
    "Encircling: By doubleckling or pressing 'd' you start drawing. Encircle the area/structure you want to annotate by simply moving
      the cursor. By double clicking again or pressing 'e' you stop drawing. Depending on the drawing mode you have chosen (Single
      or Multiple) the encircled area is highlighted immediately (Multiple) or you need to click on 'Highlight' or press 'h' (Single).",
    "",
    "Naming: After clicking on 'Highlight' you can check if the highlighted area covers the region you want to annotate.
     You are then prompted to choose the name you want to annotate the barcode spots with that fall into this area.
     This can either be a new name or one that has already been assigned within the variable. Then click on 'Name'.
     The 'Overview'-plot on the left should now display the named region in addition to all the other regions that
     you have annotated already. Naming barcode spots that have already been named results in overwriting the previous name."
  ),
  plot_orientation = create_image_annotations_descr$tab_panel_orientation,
  plot_overview = c(
    "Choose the segmentation variable that you want to alter. If you want to create a new one
    click on 'Create new segmentation variable' and you are prompted to enter the name of the new variable.
    It can then be selected. The surface plot below shows the segmentation variable you are currently working on
    and colors it according to the regions you have named.",
    "",
    "('unnamed' is the default group name that is assigned to every barcode spot by creating a new segmentation variable.)"),
  pointsize = c("The size with which points are displayed."),
  transparency_point = create_image_annotations_descr$transparency_point,
  zooming_options = c(
    "Brush the area on the plot you want to zoom in on. Then click on 'Zoom in'. You can zoom stepwise. To zoom one step
    back click on 'Zoom back'. To zoom out completely click on 'Zoom out'. Note that you can not zoom if
    you are drawing. If you want to stop drawing to zoom in on the image exit the drawing mode via shortcut 'e'
    zoom in and then start drawing again via doubleclicking or pressing 'd'. This is only possible if
    you are using drawing mode 'Single'."
  )

)

# createSpatialTrajectories()

create_spatial_trajectories_descr <- list(
  ncol = create_image_annotations_descr$ncol,
  nrow = create_image_annotations_descr$nrow,
  sgmt_size = "Size of the trajectory.",
  trajectory_ids = "Choose the IDs of the trajectories you want to plot.",
  transparency_1 = "The transparency of the points that are not included by the trajectory.",
  transparency_2 = "The transparency of the points that are included by the trajectory."
)



# all
text <- list(
  createImageAnnotations = create_image_annotations_descr,
  createSegmentation = create_segmentation_descr,
  createSpatialTrajectories = create_spatial_trajectories_descr
)
