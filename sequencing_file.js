var sequencing_file = {
    /* State needed to maintain the files table. */
    // All of the table HTML for each page in the table.
    file_table_html_by_page: {},
    // The current page number.
    current_page: null,

    document_ready: function() {
	// This function should be called as soon as the document is ready.

	// Get the csrf token so this page can post.
        $.ajaxSetup({
            data: {csrfmiddlewaretoken: Cookies.get('csrftoken') }
        });

	// Add a click event to the Files tab
	$('a[href="#files"]').on("click", function(e) {
	    sequencing_file.display_files_tab();
	    return true;
	});
    },

    display_files_tab: function() {
	/* Populate and display the files tab. */
	// If the loading pacifier exists, then the table has already been loaded.
	if ( $("#files-pacifier").length == 0) {
	    return true;
	}

	// Get the Sequence ID.
	var sequence_id = $('#files').attr('sequence_id');

	// Make an AJAX call to generate the table.
	$.get(
	    "/patient/sequence/get_sequence_file_table/" + sequence_id + "/",
	    {},
	    function (data) {
		// Display all messages to the user.
		update_messages(data.messages);

		// If the call was successful,
		if (data.success) {
		    // Remove the loading pacifier.
		    $("#files-pacifier").remove();

		    // If there is table data,
		    if (data.table_html_by_page) {
			// If there is any table data, Start on page 0
			sequencing_file.current_page = 0;
			sequencing_file.file_table_html_by_page = data.table_html_by_page;

			// Replace the div's HTML with page 0.
			$("#file_table").html(sequencing_file.file_table_html_by_page[sequencing_file.current_page]);

			// Remove the default navigation window.
			change_ability_visibility("#file_table .pagination", false, false);

			page_count = Object.keys(sequencing_file.file_table_html_by_page).length;
			// If this is the only page, hide the navigation buttons
			if (page_count <= 1) {
			    change_ability_visibility("#file_table_nav", false, false);
			}
			else {
			    // Override the navigation button behavior.
			    $("#first_page").off("click");
			    $("#first_page").click(sequencing_file.goto_first_page_of_file_table);

			    $("#previous_page").off("click");
			    $("#previous_page").click(sequencing_file.goto_previous_page_of_file_table);

			    $("#last_page").off("click");
			    $("#last_page").click(sequencing_file.goto_last_page_of_file_table);

			    $("#next_page").off("click");
			    $("#next_page").click(sequencing_file.goto_next_page_of_file_table);

			    sequencing_file.update_js_file_table_button_behavior();
			}
		    }
		    else {
			// No table data? Clear the state variables
			$("#file_table").html("No Sequencing Files found.");
		    }
		}
	    }
	);

	return true;
    },

    update_js_file_table_button_behavior: function() {
	/*
	  Update the button behavior based on the current page.
	*/

	// Gather the information needed to determine button state.
	// Are we on the last page?
	page_count = Object.keys(sequencing_file.file_table_html_by_page).length;
	on_last_page = (sequencing_file.current_page >= page_count - 1);
	// Are we on the first page?
	on_first_page = (sequencing_file.current_page == 0);

	// Previous Page button
	previous_button_css = "#file_table_nav #previous_page";
	// Show and Disable if we're on the first page.
	if (on_first_page) {
	    change_ability_visibility(previous_button_css, false, true);
	}
	else {
	    // Otherwise, Show and Enable
	    change_ability_visibility(previous_button_css, true, true);
	}

	// Next Page button
	next_button_css = "#file_table_nav #next_page";
	// Show and Disable if we're on the last page.
	if (on_last_page) {
	    change_ability_visibility(next_button_css, false, true);
	}
	else {
	    // Otherwise, Show and Enable
	    change_ability_visibility(next_button_css, true, true);
	}

	// First Page button
	first_button_css = "#file_table_nav #first_page";
	// Show and Disable if we're on the first page.
	if (on_first_page) {
	    change_ability_visibility(first_button_css, false, true);
	}
	else {
	    // Otherwise, Show and Enable
	    change_ability_visibility(first_button_css, true, true);
	}

	// Last Page button
	last_button_css = "#file_table_nav #last_page";
	// Show and Disable if we're on the last page and all of the data has been loaded.
	if (on_last_page) {
	    change_ability_visibility(last_button_css, false, true);
	}
	else {
	    // Show and Enable if we're not on the last page and all of the data has been loaded.
	    change_ability_visibility(last_button_css, true, true);
	}

	// Hide the built in navigation
	change_ability_visibility("#file_table .pagination", false, false);

	// Update the number of pages
	current_page = sequencing_file.current_page;
	displayed_page = current_page + 1;
	$('#file_table_nav #page_display').html("Page " + displayed_page + " out of " + page_count);
	change_ability_visibility("#file_table_nav #page_display", true, true);

	// Now reveal the button navigation.
	change_ability_visibility("#file_table_nav", true, true);
    },

    goto_next_page_of_file_table: function(e) {
	/*
	  Go to the next page of the file table.
	*/

	// If we're on the last page, return
	if (sequencing_file.current_page >= Object.keys(sequencing_file.file_table_html_by_page).length - 1) {
	    return;
	}

	// Update the current page number.
	sequencing_file.current_page += 1;

	// Replace the table HTML with the previous page.
	$('#file_table').html(sequencing_file.file_table_html_by_page[sequencing_file.current_page]);

	// Update the table button behavior.
	sequencing_file.update_js_file_table_button_behavior();
    },

    goto_first_page_of_file_table: function(e) {
	/*
	  Go to the first page of the gene table.
	*/

	// If we're on the first page, return
	if (sequencing_file.current_page <= 0) {
	    return;
	}

	// Update the current page number.
	sequencing_file.current_page = 0;

	// Replace the table HTML with the first page.
	$('#file_table').html(sequencing_file.file_table_html_by_page[sequencing_file.current_page]);

	// Update the table button behavior.
	sequencing_file.update_js_file_table_button_behavior();
    },

    goto_last_page_of_file_table: function(e) {
	/*
	  Go to the last page of the gene table.
	*/

	// If we're on the last page, return
	last_page_index = Object.keys(sequencing_file.file_table_html_by_page).length - 1;
	if (sequencing_file.current_page >= last_page_index) {
	    return;
	}

	// Update the current page number.
	sequencing_file.current_page = last_page_index;

	// Replace the table HTML with the previous page.
	$('#file_table').html(sequencing_file.file_table_html_by_page[sequencing_file.current_page]);

	// Update the table button behavior.
	sequencing_file.update_js_file_table_button_behavior();
    },

    goto_previous_page_of_file_table: function(e) {
	/*
	  Go to the previous page of the gene table.
	*/

	// If we're on the first page, return
	if (sequencing_file.current_page <= 0) {
	    return;
	}

	// Update the current page number.
	sequencing_file.current_page -= 1;

	// Replace the table HTML with the previous page.
	$('#file_table').html(sequencing_file.file_table_html_by_page[sequencing_file.current_page]);

	// Update the table button behavior.
	sequencing_file.update_js_file_table_button_behavior();
    },

    show_file_download_modal: function(event, element) {
	/* Populate and load a modal that will be used to download the image. */
	// Make an AJAX call to get the URL.
	var file_uuid = $(element).attr('file_uuid');
	var file_filename = $(element).attr('filename');
	var download_url = "/patient/downloadsequenceurl/" + file_uuid + "/";

	$.get(
	    download_url,
	    {},
	    function (data) {
		// Show all messages to the user.
		update_messages(data.messages);

		// If successful
		if (data.success) {
		    // Change the download url.
		    let download_url = data.download_url;
		    $('#sequence_file_download_modal .sequencing_file_url').html(download_url);
		    $('#direct_file_url').attr("href", download_url);

		    // Change the filename.
		    $('#sequence_file_download_modal .sequencing_file_name').html(file_filename);

		    // Show the modal.
		    $('#sequence_file_download_modal').modal('show');
		}
	    }
	);
    },
};

function update_messages(messages) {
    // Notify the user of all messages.
    if (typeof messages != 'undefined') {
        messages.forEach(function(entry){
            new PNotify(
                {
                    text: entry.message,
                    type: entry.extra_tags,
                    styling: 'bootstrap3'
                }
            )
        });
    }
};

function change_ability_visibility(css, is_able, is_visible) {
    // Changes whether the element identified by css is disabled and/or hidden.

    // Change ability.
    if (is_able) {
	$(css).removeClass('disabled');
    }
    else {
	$(css).addClass('disabled');
    }

    // Change visibility.
    if (is_visible) {
	$(css).removeClass('hidden');
    }
    else {
	$(css).addClass('hidden');
    }
};
