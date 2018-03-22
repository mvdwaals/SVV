

figure(666);

change_active_screen = [];
change_display_graph_state = [];

for t=2:length(flightdata.time.data)
    if flightdata.display_active_screen.data(t-1) ~= flightdata.display_active_screen.data(t)
        change_active_screen = cat(1, change_active_screen, flightdata.time.data(t));
    end
    if flightdata.display_graph_state.data(t-1) ~= flightdata.display_graph_state.data(t)
        change_display_graph_state = cat(1, change_display_graph_state, flightdata.time.data(t));
    end
end

