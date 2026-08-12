using System;
using System.Drawing;
using System.Windows.Forms;
using Rhino;
using Rhino.DocObjects;

namespace Toplet_v0_Alpha.TopOpt3D
{
    internal sealed class ResultViewerForm : Form
    {
        private readonly TopletDisplayConduit _conduit;
        private readonly RhinoDoc             _doc;
        private readonly TrackBar             _thresholdBar;
        private readonly Label                _thresholdLabel;

        public ResultViewerForm(TopletDisplayConduit conduit, RhinoDoc doc)
        {
            _conduit = conduit;
            _doc     = doc;

            Text            = "Toplet Result";
            FormBorderStyle = FormBorderStyle.FixedToolWindow;
            StartPosition   = FormStartPosition.Manual;
            Location        = new Point(60, 60);
            ClientSize      = new Size(230, 230);
            TopMost         = true;

            // Mode group
            var modeGroup = new GroupBox {
                Text     = "Display Mode",
                Location = new Point(8, 6),
                Size     = new Size(214, 108)
            };

            var rbHeatmap  = new RadioButton { Text = "Density Heatmap", Location = new Point(10, 22), AutoSize = true, Checked = true };
            var rbSolidVoid = new RadioButton { Text = "Solid / Void",    Location = new Point(10, 50), AutoSize = true };
            var rbGhosted  = new RadioButton { Text = "Ghosted",          Location = new Point(10, 78), AutoSize = true };

            rbHeatmap.CheckedChanged  += (s, e) => { if (rbHeatmap.Checked)   SetMode(TopletDisplayMode.DensityHeatmap); };
            rbSolidVoid.CheckedChanged += (s, e) => { if (rbSolidVoid.Checked) SetMode(TopletDisplayMode.SolidVoid); };
            rbGhosted.CheckedChanged  += (s, e) => { if (rbGhosted.Checked)   SetMode(TopletDisplayMode.Ghosted); };

            modeGroup.Controls.AddRange(new Control[] { rbHeatmap, rbSolidVoid, rbGhosted });
            Controls.Add(modeGroup);

            // Threshold row
            Controls.Add(new Label { Text = "Threshold:", Location = new Point(10, 122), AutoSize = true });
            _thresholdLabel = new Label { Text = "0.50", Location = new Point(190, 122), AutoSize = true };
            Controls.Add(_thresholdLabel);

            _thresholdBar = new TrackBar {
                Minimum       = 1,
                Maximum       = 99,
                Value         = 50,
                TickFrequency = 10,
                Location      = new Point(8, 138),
                Size          = new Size(214, 32)
            };
            _thresholdBar.ValueChanged += OnThresholdChanged;
            Controls.Add(_thresholdBar);

            // Buttons
            var btnAdd = new Button { Text = "Add to Document", Location = new Point(8,  188), Size = new Size(130, 28) };
            var btnClose = new Button { Text = "Close",          Location = new Point(146, 188), Size = new Size(76,  28) };

            btnAdd.Click   += OnAddToDocument;
            btnClose.Click += (s, e) => Close();
            Controls.AddRange(new Control[] { btnAdd, btnClose });

            FormClosing += (s, e) => {
                _conduit.Enabled = false;
                _doc.Views.Redraw();
            };
        }

        private void SetMode(TopletDisplayMode mode)
        {
            _conduit.Mode = mode;
            _doc.Views.Redraw();
        }

        private void OnThresholdChanged(object sender, EventArgs e)
        {
            double t = _thresholdBar.Value / 100.0;
            _thresholdLabel.Text = $"{t:F2}";
            _conduit.Threshold   = t;
            _doc.Views.Redraw();
        }

        private void OnAddToDocument(object sender, EventArgs e)
        {
            Rhino.Geometry.Mesh mesh = _conduit.BakeMesh();
            if (mesh == null || !mesh.IsValid)
            {
                MessageBox.Show("Nothing to add — try adjusting the threshold.", "Toplet");
                return;
            }

            int layerIdx = Visualization3D.EnsureLayer(_doc, "Toplet3D_Result");
            var attr = new ObjectAttributes { LayerIndex = layerIdx };
            _doc.Objects.AddMesh(mesh, attr);
            _doc.Views.Redraw();
        }
    }
}
