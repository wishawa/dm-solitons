from plotly.subplots import make_subplots
import plotly.graph_objects as pgo
I3D = {"is_3d": True}
def plot_multiple_3Ds(plots):
	rows = len(plots)
	cols = max([len(ro) for ro in plots])
	fig = make_subplots(
		rows = rows,
		cols = cols,
		specs=[[I3D] * cols] * rows
	)
	# fig.update_layout(height = 500, width = 1500)
	for (i,ro) in enumerate(plots):
		for (j,pl) in enumerate(ro):
			fig.add_trace(pl, i+1, j+1)
	fig = pgo.Figure(fig)

	js = '''
	(function() {
		const gd = document.getElementById('{plot_id}');
		let isUnderRelayout = 0;

		gd.on('plotly_relayout', (e) => {
			const camera = Object.entries(e).find(([k, v]) => k.startsWith("scene") && k.endsWith("camera"))[1];
			if (camera && isUnderRelayout == 0) {
				const relayoutKeys = Object.keys(gd.layout).filter(k => k.startsWith("scene"));
				isUnderRelayout = relayoutKeys.length;
				relayoutKeys.forEach(k => {
					Plotly.relayout(gd, k + ".camera", camera)
					.then(() => {
						isUnderRelayout--;
					});
				});
			}
		})
	})();
	'''
	return fig, js