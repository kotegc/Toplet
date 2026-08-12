using System;
using System.Drawing;
using System.Windows.Forms;

namespace Toplet_v0_Alpha.TopOpt3D
{
    internal sealed class TopletSetupForm : Form
    {
        public double CellSize      { get; private set; }
        public bool   LoadIsFace    { get; private set; }
        public bool   SupportIsFace { get; private set; }

        private readonly double  _width, _depth, _height;
        private readonly TextBox _cellSizeBox;
        private readonly Label   _previewLabel;
        private readonly Button  _okButton;

        private readonly RadioButton _rbLoadFace;
        private readonly RadioButton _rbSupportFaces;

        public TopletSetupForm(double width, double depth, double height, double defaultCellSize)
        {
            _width  = width;
            _depth  = depth;
            _height = height;

            Text            = "Toplet3D Setup";
            FormBorderStyle = FormBorderStyle.FixedDialog;
            StartPosition   = FormStartPosition.CenterScreen;
            ClientSize      = new Size(290, 330);
            MaximizeBox     = false;
            MinimizeBox     = false;
            TopMost         = true;

            // --- Element size ---
            Controls.Add(new Label {
                Text = "Element Size:",
                Location = new Point(12, 14),
                AutoSize = true
            });

            _cellSizeBox = new TextBox {
                Text     = defaultCellSize.ToString("G4"),
                Location = new Point(12, 34),
                Size     = new Size(110, 22)
            };
            _cellSizeBox.TextChanged += (s, e) => {
                if (double.TryParse(_cellSizeBox.Text, out double cs) && cs > 0)
                    UpdatePreview(cs);
            };
            Controls.Add(_cellSizeBox);

            _previewLabel = new Label {
                Location  = new Point(12, 62),
                Size      = new Size(266, 52),
                ForeColor = Color.DimGray
            };
            Controls.Add(_previewLabel);

            // --- Load type ---
            var loadGroup = new GroupBox {
                Text     = "Load Type",
                Location = new Point(12, 122),
                Size     = new Size(266, 72)
            };
            var rbLoadPoint = new RadioButton {
                Text    = "Point Load",
                Location = new Point(10, 22),
                AutoSize = true,
                Checked  = true
            };
            _rbLoadFace = new RadioButton {
                Text     = "Face Load",
                Location = new Point(10, 46),
                AutoSize = true
            };
            loadGroup.Controls.AddRange(new Control[] { rbLoadPoint, _rbLoadFace });
            Controls.Add(loadGroup);

            // --- Support type ---
            var supportGroup = new GroupBox {
                Text     = "Support Type",
                Location = new Point(12, 204),
                Size     = new Size(266, 72)
            };
            var rbSupportPoints = new RadioButton {
                Text     = "Support Points",
                Location = new Point(10, 22),
                AutoSize = true,
                Checked  = true
            };
            _rbSupportFaces = new RadioButton {
                Text     = "Support Faces",
                Location = new Point(10, 46),
                AutoSize = true
            };
            supportGroup.Controls.AddRange(new Control[] { rbSupportPoints, _rbSupportFaces });
            Controls.Add(supportGroup);

            // --- Buttons ---
            _okButton = new Button {
                Text         = "Continue",
                Location     = new Point(158, 292),
                Size         = new Size(80, 28),
                DialogResult = DialogResult.OK
            };
            var cancelButton = new Button {
                Text         = "Cancel",
                Location     = new Point(68, 292),
                Size         = new Size(80, 28),
                DialogResult = DialogResult.Cancel
            };
            Controls.AddRange(new Control[] { _okButton, cancelButton });

            // Do NOT set AcceptButton / CancelButton here — any Enter or Esc keypress
            // queued from the preceding Rhino GetOneObject call would fire them
            // immediately, closing the dialog before the user sees it.
            // Wire them up in Shown, after the message queue has been flushed.
            Shown += (s, e) => {
                AcceptButton = _okButton;
                CancelButton = cancelButton;
                Activate();
            };

            FormClosing += (s, e) => {
                if (DialogResult == DialogResult.OK)
                {
                    CellSize      = double.TryParse(_cellSizeBox.Text, out double cs) && cs > 0
                                    ? cs : defaultCellSize;
                    LoadIsFace    = _rbLoadFace.Checked;
                    SupportIsFace = _rbSupportFaces.Checked;
                }
            };

            UpdatePreview(defaultCellSize);
        }

        private void UpdatePreview(double cs)
        {
            int nelx = Math.Max(2, (int)Math.Ceiling(_width  / cs)); if (nelx % 2 != 0) nelx++;
            int nely = Math.Max(2, (int)Math.Ceiling(_depth  / cs)); if (nely % 2 != 0) nely++;
            int nelz = Math.Max(2, (int)Math.Ceiling(_height / cs)); if (nelz % 2 != 0) nelz++;

            long elems = (long)nelx * nely * nelz;
            long dofs  = (long)(nelx + 1) * (nely + 1) * (nelz + 1) * 3;

            string note = dofs > 3_000_000 ? "  (!) very large — may be slow"
                        : dofs > 1_000_000 ? "  — large grid"
                        : "";

            _previewLabel.Text = string.Format(
                "Grid:  {0} x {1} x {2}\nElements:  {3:N0}    DOFs:  {4:N0}{5}",
                nelx, nely, nelz, elems, dofs, note);

            _previewLabel.ForeColor = dofs > 3_000_000 ? Color.Firebrick
                                    : dofs > 1_000_000 ? Color.DarkOrange
                                    : Color.DimGray;

            _okButton.Enabled = dofs <= 10_000_000;
        }
    }
}
