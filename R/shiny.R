.shiny_primary_plot = function(plot_ct2, tx_df, plotType, selAttribute, selIds, plotBrush){
    #visible binding NOTE
    id = NULL
    plot_ct2 = setRegionMetaData(plot_ct2, tx_df)
    theme_update(panel.background = element_blank(), panel.grid = element_blank())
    switch(plotType,
           signal_max = {
               p = plotDimReducePoints(plot_ct2, point_size = .8)
           },
           annotation = {
               p = plotDimReducePoints(plot_ct2, color_VAR = selAttribute, point_size = .8)
           },
           signal_max_vs_annotation = {
               p = plotDimReducePoints(plot_ct2, extra_VARS = selAttribute, point_size = .8) +
                   facet_grid(paste0(getNameVariable(plot_ct2), "~", selAttribute))
           }, {
               stop("bad plotType: ", plotType)
           })
    p = p + guides(color = guide_legend(override.aes = list(size = 2.5)))
    if(!is.null(plotBrush)){
        p = p + annotate("rect",
                         xmin = plotBrush$xmin,
                         xmax = plotBrush$xmax,
                         ymin = plotBrush$ymin,
                         ymax = plotBrush$ymax,
                         fill = "pink",
                         color = "pink",
                         alpha = .2)
    }
    if(!is.null(selIds)){
        tx_df.sel = dplyr::filter(tx_df, id %in% selIds)
        p = p + annotate("point",
                         size = 2.5,
                         shape = 1,
                         x = tx_df.sel$tx,
                         y = tx_df.sel$ty,
                         color = "red")
    }

    p + coord_fixed()

}

#' .shiny_secondary_plot
#'
#' @examples
#'
#' plot_ct2 = exampleChIPtsne2.with_meta()
#' plot_ct2 = dimReduceTSNE(plot_ct2)
#' tx_df = getRegionMetaData(plot_ct2)
#' plotType = "line_color_selection"
#' selAttribute = "peak_MCF10A_CTCF"
#' selIds = character()
#' selIds = rownames(plot_ct2)[1:10]
#'
#' chiptsne2:::.shiny_secondary_plot(
#'   plot_ct2,
#'   tx_df,
#'   "line_color_selection",
#'   selAttribute = "peak_MCF10A_CTCF",
#'   selIds = rownames(plot_ct2)[1:10]
#' )
#' chiptsne2:::.shiny_secondary_plot(
#'   plot_ct2,
#'   tx_df,
#'   "line_color_selection",
#'   selAttribute = "peak_MCF10A_CTCF",
#'   selIds = character()
#' )
#' chiptsne2:::.shiny_secondary_plot(
#'   plot_ct2,
#'   tx_df,
#'   "line_color_annotation",
#'   selAttribute = "peak_MCF10A_CTCF",
#'   selIds = rownames(plot_ct2)[1:10]
#' )
#' chiptsne2:::.shiny_secondary_plot(
#'   plot_ct2,
#'   tx_df,
#'   "line_color_annotation",
#'   selAttribute = "peak_MCF10A_CTCF",
#'   selIds = character()
#' )
#' chiptsne2:::.shiny_secondary_plot(
#'   plot_ct2,
#'   tx_df,
#'   "heatmap_cluster_selection",
#'   selAttribute = "peak_MCF10A_CTCF",
#'   selIds = rownames(plot_ct2)[1:10]
#' )
#' chiptsne2:::.shiny_secondary_plot(
#'   plot_ct2,
#'   tx_df,
#'   "heatmap_cluster_selection",
#'   selAttribute = "peak_MCF10A_CTCF",
#'   selIds = character()
#' )
#' chiptsne2:::.shiny_secondary_plot(
#'   plot_ct2,
#'   tx_df,
#'   "heatmap_cluster_annotation",
#'   selAttribute = "peak_MCF10A_CTCF",
#'   selIds = rownames(plot_ct2)[1:10]
#' )
#' chiptsne2:::.shiny_secondary_plot(
#'   plot_ct2,
#'   tx_df,
#'   "heatmap_cluster_annotation",
#'   selAttribute = "peak_MCF10A_CTCF",
#'   selIds = character()
#' )
#'
.shiny_secondary_plot = function(plot_ct2, tx_df, plotType, selAttribute, selIds){
    theme_update(panel.background = element_blank(), panel.grid = element_blank())
    if(length(selIds) > 0){
        plot_FUN = .shiny_secondary_plot.has_selection

    }else{
        plot_FUN = .shiny_secondary_plot.no_selection
    }
    do.call(plot_FUN, args = get_args(to_ignore = c("plot_FUN")))
}

.shiny_secondary_plot.no_selection = function(plot_ct2, tx_df, plotType, selAttribute, selIds){
    plot_ct2 = setRegionMetaData(plot_ct2, tx_df)
    switch(plotType,
           line_color_selection = {
               p = plotSignalLinePlot(plot_ct2, group_VAR = selAttribute, n_splines = 5, moving_average_window = 1)
           },
           line_color_annotation = {
               line_colors = seqsetvis::safeBrew(as.character(tx_df[[selAttribute]]))
               p = plotSignalLinePlot(plot_ct2, color_VAR = selAttribute, n_splines = 5, moving_average_window = 1) +
                   scale_color_manual(values = line_colors)
           },
           heatmap_cluster_selection = {
               p = plotSignalHeatmap(plot_ct2, group_VARS = c(selAttribute), relative_heatmap_width = .83, relative_heatmap_height = .65) +
                   labs(subtitle = "no selection")#, annotation_colors = list(line_colors))
           },
           heatmap_cluster_annotation = {
               line_colors = seqsetvis::safeBrew(as.character(tx_df[[selAttribute]]))
               p = plotSignalHeatmap(plot_ct2, group_VARS = c(selAttribute), relative_heatmap_width = .83, relative_heatmap_height = .65)
           }, {
               stop("bad plotType: ", plotType)
           }
    )
    p
}
.shiny_secondary_plot.has_selection = function(plot_ct2, tx_df, plotType, selAttribute, selIds){
    tx_df$is_selected = "not selected"
    tx_df$is_selected[tx_df[[plot_ct2@region_VAR]] %in% selIds] = "selected"
    plot_ct2 = setRegionMetaData(plot_ct2, tx_df)
    switch(plotType,
           line_color_selection = {
               p = plotSignalLinePlot(plot_ct2, group_VAR = selAttribute, color_VAR = "is_selected", n_splines = 5, moving_average_window = 1) +
                   scale_color_manual(values = c("not selected" = 'gray60', "selected" = "pink"))
           },
           line_color_annotation = {
               line_colors = seqsetvis::safeBrew(as.character(tx_df[[selAttribute]]))
               p = plotSignalLinePlot(plot_ct2, color_VAR = selAttribute, group_VAR = "is_selected", n_splines = 5, moving_average_window = 1) +
                   scale_color_manual(values = line_colors)
           },
           heatmap_cluster_selection = {
               p = plotSignalHeatmap(plot_ct2, group_VARS = c(selAttribute, "is_selected"), relative_heatmap_width = .66, relative_heatmap_height = .65)
           },
           heatmap_cluster_annotation = {
               p = plotSignalHeatmap(plot_ct2, group_VARS = c("is_selected", selAttribute), relative_heatmap_width = .66, relative_heatmap_height = .65)
           }, {
               stop("bad plotType: ", plotType)
           }
    )
    p
}

#' run_ChIPtsne2_shiny
#'
#' @param ct2 `r doc_ct2_nrr()`
#'
#' @return ChIPtsne2 object modified with user's interactive annotations.
#' @export
#' @import shiny
#'
#' @examples
#' ct2 = exampleChIPtsne2.with_meta()
#' ct2 = dimReduceTSNE(ct2)
#' \dontrun{
#'   plotDimReducePoints(ct2, color_VAR = "peak_MCF10A_CTCF")
#'   ct2.ann = run_ChIPtsne2_shiny(ct2)
#'   plotDimReducePoints(ct2.ann, color_VAR = "peak_MCF10A_CTCF")
#' }
#'
#'
run_ChIPtsne2_shiny = function(ct2){
    #visible binding NOTE
    dist_ = rnk = tx = ty = vals = NULL
    theme_update(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
    )
    # theme(axis.title = element_text(size = 14), legend.text = )
    meta_df = getRegionMetaData(ct2)
    cn = colnames(meta_df)
    num_cn = sapply(meta_df, is.numeric)
    ko_cn = c("tx", "ty", ct2@region_VAR, colnames(meta_df)[num_cn])
    cn = setdiff(cn, ko_cn)
    shiny_output = NULL

    ui = fluidPage(
        tags$h3("ChIPtsne2"),
        tags$figcaption('click "finish" to return a ChIPtsne2 with your modifications'),
        actionButton("btn_finish", label = "Finish"),
        tags$figcaption("add new annotation attribute"),
        actionButton("btn_addAttribute", label = "New Attribute"),
        tags$figcaption("select annotation attribute and add new value categories"),
        selectizeInput("txt_selAttribute", choices = cn, label = "Select Attribute"),
        textInput("txt_newValue", label = "New Category Value", value = "", placeholder = "annotation value"),
        actionButton("btn_addAnnotation", label = "Annotate Selected"),
        tags$figcaption("double-click will select the nearest n points"),
        numericInput("num_nClosest", label = "Neareset Points", min = 1, max = nrow(meta_df), value = round(.1*nrow(meta_df))),
        tags$figcaption("select points below to annotate. brush or double-click to select"),
        plotOutput(
            'plot',
            click = "plot_click",
            dblclick = "plot_dblclick",
            hover = "plot_hover",
            brush = "plot_brush"
        ),
        radioButtons(
            "txt_plotType",
            label = "Primary Plot Type",
            choiceNames = c("signal max", "annotation", "signal max / annotation"),
            choiceValues = c("signal_max", "annotation", "signal_max_vs_annotation")
        ),
        plotOutput('plot2'),
        radioButtons(
            "txt_plot2Type",
            label = "Secondary Plot Type",
            choiceNames = c("line plot selection", "line plot annotation", "heatmap selection / annotation", "heatmap annotation / selection"),
            choiceValues = c("line_color_selection", "line_color_annotation", "heatmap_cluster_selection", "heatmap_cluster_annotation")
        ),
        verbatimTextOutput("info")
    )

    server = function(input, output){
        output$plot <- renderPlot({
            .shiny_primary_plot(plot_ct2 = ct2, tx_df = r_meta_df(), input$txt_plotType, input$txt_selAttribute,  r_sel_ids(), isolate(input$plot_brush))
        })
        output$plot2 = renderPlot({
            .shiny_secondary_plot(plot_ct2 = ct2, tx_df = r_meta_df(), input$txt_plot2Type, input$txt_selAttribute,  r_sel_ids())
        })

        output$info <- renderText({
            xy_str <- function(e) {
                if(is.null(e)) return("NULL")
                paste0("x=", round(e$x, 1), " y=", round(e$y, 1))
            }
            xy_range_str <- function(e) {
                if(is.null(e)) return("NULL")
                paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1),
                       " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
            }

            m_df = r_meta_df()
            tab_att = table(m_df[, input$txt_selAttribute])
            att_str = paste(paste0(names(tab_att), ":", tab_att), collapse = ", ")
            paste(
                sep = "\n",
                paste0("click: ", xy_str(input$plot_click)),
                paste0("dblclick: ", xy_str(input$plot_dblclick)),
                paste0("hover: ", xy_str(input$plot_hover)),
                paste0("brush: ", xy_range_str(input$plot_brush)),
                paste0("selected: ", paste(r_sel_ids(), collapse = ", ")),
                paste0("annotation: ", att_str)
            )
        })

        r_sel_ids = reactiveVal(NULL)
        r_meta_df = reactiveVal(meta_df)
        r_cn = reactiveVal(colnames(meta_df))
        r_new_cn = reactiveVal(NULL)

        observeEvent({
            input$plot_dblclick
        }, {
            np_df = nearPoints(df = meta_df,
                               input$plot_dblclick,
                               xvar = "tx",
                               yvar = "ty",
                               addDist = TRUE,
                               maxpoints = 10,
                               allRows = TRUE)
            np_df.sel = np_df %>%
                dplyr::mutate(rnk = rank(dist_, ties.method = "first")) %>%
                dplyr::filter(rnk <= input$num_nClosest)
            r_sel_ids(np_df.sel$id)

        })

        observeEvent({
            input$plot_brush
        }, {
            bp_df.sel = meta_df %>% dplyr::filter(tx >= input$plot_brush$xmin & tx <= input$plot_brush$xmax &
                                                      ty >= input$plot_brush$ymin & ty <= input$plot_brush$ymax)
            r_sel_ids(bp_df.sel$id)
        })

        # Return the UI for a modal dialog with data selection input. If 'failed' is
        # TRUE, then display a message that the previous value was invalid.
        dataModal <- function(failed = FALSE) {
            modalDialog(
                textInput("newAttribute", "Attribute Name",
                          placeholder = 'name for the new attribute.'
                ),
                textInput("defaultValue", "Default Value",
                          placeholder = 'starting value for this new attribute.'
                ),
                span('New attribute has to be not present already.'),
                if (failed)
                    div(tags$b("Invalid name for new attribute.", style = "color: red;")),

                footer = tagList(
                    modalButton("Cancel"),
                    actionButton("ok", "OK")
                )
            )
        }

        observeEvent({
            input$btn_addAttribute
        }, {
            showModal(dataModal())
        })

        # When OK button is pressed, attempt to load the data set. If successful,
        # remove the modal. If not show another modal, but this time with a failure
        # message.
        observeEvent(input$ok, {
            # Check that data object exists and is data frame.
            if (!input$newAttribute %in% colnames(meta_df)){
                new_name = input$newAttribute
                def_value = input$defaultValue
                m_df = r_meta_df()
                m_df[[new_name]] = def_value
                r_meta_df(m_df)
                removeModal()
            } else {
                showModal(dataModal(failed = TRUE))
            }
        })

        # Display information about selected data
        output$dataInfo <- renderPrint({
            if (is.null(vals$data))
                "No data selected"
            else
                summary(vals$data)
        })

        observeEvent({
            r_cn() #only update when colnames actually change
        },{
            updateSelectInput(inputId = "txt_selAttribute", choices = setdiff(colnames(r_meta_df()), ko_cn), selected = r_new_cn())
        })

        observeEvent({
            r_meta_df()
        }, {
            old_cn = r_cn()
            new_cn = colnames(r_meta_df())
            if(!setequal(old_cn, new_cn)){
                r_cn(new_cn)
                r_new_cn(setdiff(new_cn, old_cn))
            }
        })

        observeEvent({
            input$btn_addAnnotation
        }, {
            m_df = r_meta_df()
            apply_annotation = TRUE
            if(input$txt_newValue == ""){
                showNotification("No Annotation Value specified, not effect.", type = "warning")
                apply_annotation = FALSE
            }
            if(length(r_sel_ids()) == 0){
                showNotification("No regions selected for annotation, no effect.", type = "warning")
                apply_annotation = FALSE
            }
            if(apply_annotation){
                m_df[r_sel_ids(), input$txt_selAttribute] = input$txt_newValue
                r_meta_df(m_df)
            }
        })

        observeEvent({
            input$btn_finish
        }, {
            new_meta = r_meta_df()
            new_ct2 = setRegionMetaData(ct2, new_meta)
            stopApp(returnValue = new_ct2)
        })

    }
    app = shinyApp(
        ui,
        server
    )
    runApp(app)
}

