using System;
using System.Drawing;
using System.Windows.Forms;

namespace Toplet_v0_Alpha.Interop
{
    internal sealed class SolverProgressForm : Form
    {
        private readonly ProgressBar _bar;
        private readonly Label _iterLabel;
        private readonly Label _compLabel;

        public SolverProgressForm(int maxIter)
        {
            Text = "Topology Optimization";
            FormBorderStyle = FormBorderStyle.FixedDialog;
            MaximizeBox = false;
            MinimizeBox = false;
            ControlBox = false;
            StartPosition = FormStartPosition.CenterScreen;
            ClientSize = new Size(380, 110);

            Controls.Add(new Label {
                Text = "Solving, please wait...",
                AutoSize = true,
                Location = new Point(12, 12),
                Font = new Font(Font, FontStyle.Bold)
            });

            _bar = new ProgressBar {
                Minimum = 0,
                Maximum = Math.Max(maxIter, 1),
                Value = 0,
                Location = new Point(12, 38),
                Size = new Size(356, 20),
                Style = ProgressBarStyle.Continuous
            };
            Controls.Add(_bar);

            _iterLabel = new Label {
                Text = $"Iteration: 0 / {maxIter}",
                AutoSize = true,
                Location = new Point(12, 70)
            };
            Controls.Add(_iterLabel);

            _compLabel = new Label {
                Text = "Compliance: —",
                AutoSize = true,
                Location = new Point(210, 70)
            };
            Controls.Add(_compLabel);
        }

        public void UpdateProgress(int iter, int maxIter, double compliance)
        {
            _bar.Value = Math.Min(iter, _bar.Maximum);
            _iterLabel.Text = $"Iteration: {iter} / {maxIter}";
            _compLabel.Text = $"Compliance: {compliance:F4}";
        }
    }

    internal sealed class SolverCompletedForm : Form
    {
        public SolverCompletedForm(int iterations, int maxIter, double compliance, TimeSpan elapsed, bool converged)
        {
            Text = "Optimization Complete";
            FormBorderStyle = FormBorderStyle.FixedDialog;
            MaximizeBox = false;
            MinimizeBox = false;
            StartPosition = FormStartPosition.CenterScreen;
            ClientSize = new Size(310, 158);

            string status = converged ? "Converged" : "Max iterations reached";
            var rows = new[]
            {
                $"Status:       {status}",
                $"Iterations:   {iterations} / {maxIter}",
                $"Compliance:   {compliance:F6}",
                $"Elapsed:      {elapsed.TotalSeconds:F2} s",
            };

            var mono = new Font("Consolas", 9.5f);
            int y = 14;
            foreach (var row in rows)
            {
                Controls.Add(new Label { Text = row, AutoSize = true, Location = new Point(16, y), Font = mono });
                y += 24;
            }

            var ok = new Button {
                Text = "OK",
                DialogResult = DialogResult.OK,
                Location = new Point(214, y + 2),
                Size = new Size(80, 26)
            };
            Controls.Add(ok);
            AcceptButton = ok;
        }
    }
}
