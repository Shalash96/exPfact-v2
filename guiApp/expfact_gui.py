import sys
import os
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QTabWidget, QWidget, QVBoxLayout, QHBoxLayout,
    QFormLayout, QLabel, QLineEdit, QPushButton, QFileDialog, QSpinBox,
    QDoubleSpinBox, QCheckBox, QComboBox, QTextEdit, QMessageBox, QStatusBar
)
from PyQt6.QtCore import QProcess, Qt
from PyQt6.QtGui import QPixmap

# --- Configuration ---
GUI_APP_LOCATION = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.join(os.path.dirname(GUI_APP_LOCATION), 'python')
R_SCRIPT_DIR_RELATIVE_TO_PY_SCRIPT = "../R" # Relative to Python scripts in SCRIPT_DIR

# --- Change CWD to project root, then to 'results' subdirectory ---
PROJECT_ROOT = os.path.dirname(GUI_APP_LOCATION)
os.chdir(PROJECT_ROOT)
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)
os.chdir(RESULTS_DIR) # Scripts will run with 'results' as CWD

# --- Helper function to get full script path ---
def get_script_path(script_name):
    return os.path.join(SCRIPT_DIR, script_name)

class ExPfactGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ExPfact Suite GUI")
        self.setGeometry(100, 100, 900, 700)

        self.process = None

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QVBoxLayout(self.central_widget)

        self.tab_widget = QTabWidget()
        self.main_layout.addWidget(self.tab_widget)

        self._create_expfact_tab()
        self._create_md2pfact_tab()
        self._create_pfact2dpred_tab()
        self._create_isotopes_tab()
        self._create_utilities_tab()

        self.log_output_area = QTextEdit()
        self.log_output_area.setReadOnly(True)
        self.main_layout.addWidget(self.log_output_area, 1)

        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("ExPfact Suite GUI Ready. Devloped by: Mahmoud Shalash a member of Paci Lab, University of Bologna, Italy.")


    def _log_message(self, message, is_error=False):
        if is_error:
            self.log_output_area.append(f"<font color='red'>{message.strip()}</font>")
        else:
            self.log_output_area.append(message.strip())
        self.log_output_area.ensureCursorVisible()

    def _browse_file(self, line_edit, caption="Open File", file_filter="All Files (*.*)"):
        file_name, _ = QFileDialog.getOpenFileName(self, caption, "", file_filter)
        if file_name:
            line_edit.setText(file_name)

    def _browse_directory(self, line_edit, caption="Select Directory"):
        dir_name = QFileDialog.getExistingDirectory(self, caption)
        if dir_name:
            line_edit.setText(dir_name)

    def _run_script(self, script_name, args_list, script_is_python=True):
        if self.process and self.process.state() == QProcess.ProcessState.Running:
            QMessageBox.warning(self, "Process Busy", "Another process is already running.")
            return

        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._handle_stdout)
        self.process.finished.connect(self._process_finished)

        executable = "python"
        full_script_path_or_command = get_script_path(script_name)

        if script_is_python and not os.path.exists(full_script_path_or_command):
            self._log_message(f"ERROR: Python script not found: {full_script_path_or_command}", is_error=True)
            self.status_bar.showMessage(f"Error: Script not found: {script_name}")
            self.process = None
            return
        
        command_display_parts = []
        effective_args = []

        if script_is_python:
            executable = "python"
            command_display_parts.append(executable)
            command_display_parts.append(full_script_path_or_command)
            command_display_parts.extend(args_list)
            effective_args = [full_script_path_or_command] + args_list
        else: 
            executable = script_name # e.g., "Rscript"
            command_display_parts.append(executable)
            command_display_parts.extend(args_list)
            effective_args = args_list


        self._log_message(f"Running: {' '.join(command_display_parts)}")
        self.status_bar.showMessage(f"Running {script_name if script_is_python else executable}...")

        self._toggle_run_buttons(enable=False)
        self.process.start(executable, effective_args)


    def _toggle_run_buttons(self, enable=True):
        for i in range(self.tab_widget.count()):
            tab_content = self.tab_widget.widget(i)
            sub_tabs = tab_content.findChildren(QTabWidget)
            if sub_tabs:
                for sub_tab_widget in sub_tabs:
                    for j in range(sub_tab_widget.count()):
                        sub_tab_content = sub_tab_widget.widget(j)
                        buttons = sub_tab_content.findChildren(QPushButton)
                        for button in buttons:
                            if "run_" in button.objectName().lower() or \
                               "calculate" in button.text().lower() or \
                               "process" in button.text().lower():
                                button.setEnabled(enable)
            else:
                buttons = tab_content.findChildren(QPushButton)
                for button in buttons:
                    if "run_" in button.objectName().lower() or \
                       "calculate" in button.text().lower() or \
                       "process" in button.text().lower():
                        button.setEnabled(enable)

    def _handle_stdout(self):
        data = self.process.readAllStandardOutput().data().decode().strip()
        if data:
            self._log_message(data)

    def _process_finished(self, exit_code, exit_status):
        self._toggle_run_buttons(enable=True)

        if exit_status == QProcess.ExitStatus.CrashExit:
            self._log_message("Process crashed.", is_error=True)
            self.status_bar.showMessage("Process crashed.")
        elif exit_code != 0:
            self._log_message(f"Process finished with error code: {exit_code}.", is_error=True)
            self.status_bar.showMessage(f"Process failed (code: {exit_code}).")
        else:
            self._log_message("Process finished successfully.")
            self.status_bar.showMessage("Process completed.")

        active_tab_text = self.tab_widget.tabText(self.tab_widget.currentIndex())
        if active_tab_text == "MD to Pfact" and hasattr(self, 'md2pfact_out_edit') and \
           self.md2pfact_out_edit.text() and exit_code == 0:
            plot_path = f"{self.md2pfact_out_edit.text()}_avg.png"
            if os.path.exists(plot_path):
                self._display_md_plot(plot_path)
            else:
                self._log_message(f"Plot not found: {plot_path}", is_error=True)
        self.process = None

    def _create_expfact_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        form_layout = QFormLayout()

        self.expfact_dexp_edit = QLineEdit()
        self.expfact_dexp_browse = QPushButton("Browse...")
        self.expfact_dexp_browse.clicked.connect(lambda: self._browse_file(self.expfact_dexp_edit, "Dexp File", "Data Files (*.dexp *.txt *.dat);;All Files (*.*)"))
        hb_dexp = QHBoxLayout()
        hb_dexp.addWidget(self.expfact_dexp_edit)
        hb_dexp.addWidget(self.expfact_dexp_browse)
        form_layout.addRow("Dexp File (--dexp):", hb_dexp)

        self.expfact_ass_edit = QLineEdit()
        self.expfact_ass_browse = QPushButton("Browse...")
        self.expfact_ass_browse.clicked.connect(lambda: self._browse_file(self.expfact_ass_edit, "Assignments File", "List Files (*.ass *.list *.txt);;All Files (*.*)"))
        hb_ass = QHBoxLayout()
        hb_ass.addWidget(self.expfact_ass_edit)
        hb_ass.addWidget(self.expfact_ass_browse)
        form_layout.addRow("Assignments File (--ass):", hb_ass)

        self.expfact_seq_edit = QLineEdit()
        self.expfact_seq_browse = QPushButton("Browse...")
        self.expfact_seq_browse.clicked.connect(lambda: self._browse_file(self.expfact_seq_edit, "Sequence File", "Sequence Files (*.seq *.fasta *.txt);;All Files (*.*)"))
        hb_seq = QHBoxLayout()
        hb_seq.addWidget(self.expfact_seq_edit)
        hb_seq.addWidget(self.expfact_seq_browse)
        form_layout.addRow("Sequence File (--seq):", hb_seq)

        self.expfact_pfact_cb = QCheckBox("Enable Initial Pfact File")
        self.expfact_pfact_edit = QLineEdit()
        self.expfact_pfact_edit.setEnabled(False)
        self.expfact_pfact_browse = QPushButton("Browse...")
        self.expfact_pfact_browse.setEnabled(False)
        self.expfact_pfact_cb.toggled.connect(self.expfact_pfact_edit.setEnabled)
        self.expfact_pfact_cb.toggled.connect(self.expfact_pfact_browse.setEnabled)
        self.expfact_pfact_browse.clicked.connect(lambda: self._browse_file(self.expfact_pfact_edit, "Pfact File", "Pfact Files (*.pfact *.txt);;All Files (*.*)"))
        hb_pfact = QHBoxLayout()
        hb_pfact.addWidget(self.expfact_pfact_edit)
        hb_pfact.addWidget(self.expfact_pfact_browse)
        form_layout.addRow(self.expfact_pfact_cb, hb_pfact)

        self.expfact_temp_spin = QDoubleSpinBox()
        self.expfact_temp_spin.setRange(0, 500)
        self.expfact_temp_spin.setValue(300)
        form_layout.addRow("Temperature (--temp, K):", self.expfact_temp_spin)

        self.expfact_ph_spin = QDoubleSpinBox()
        self.expfact_ph_spin.setRange(0, 14)
        self.expfact_ph_spin.setValue(7.0)
        form_layout.addRow("pH (--pH):", self.expfact_ph_spin)

        self.expfact_penalty_spin = QDoubleSpinBox()
        self.expfact_penalty_spin.setRange(0, 100)
        self.expfact_penalty_spin.setValue(0.0)
        form_layout.addRow("Penalty Factor (--harm):", self.expfact_penalty_spin)

        self.expfact_tol_spin = QDoubleSpinBox()
        self.expfact_tol_spin.setDecimals(8)
        self.expfact_tol_spin.setSingleStep(1e-7)
        self.expfact_tol_spin.setValue(1e-6)
        form_layout.addRow("Tolerance (--tol):", self.expfact_tol_spin)

        self.expfact_rand_cb = QCheckBox("Enable Random Search Steps")
        self.expfact_rand_spin = QSpinBox()
        self.expfact_rand_spin.setRange(0, 100000)
        self.expfact_rand_spin.setValue(1000)
        self.expfact_rand_spin.setEnabled(False)
        self.expfact_rand_cb.toggled.connect(self.expfact_rand_spin.setEnabled)
        form_layout.addRow(self.expfact_rand_cb, self.expfact_rand_spin)

        self.expfact_rep_spin = QSpinBox()
        self.expfact_rep_spin.setRange(1, 100)
        self.expfact_rep_spin.setValue(1)
        form_layout.addRow("Replicates (--rep):", self.expfact_rep_spin)

        self.expfact_ncores_spin = QSpinBox()
        self.expfact_ncores_spin.setRange(1, os.cpu_count() or 1)
        self.expfact_ncores_spin.setValue(1)
        form_layout.addRow("Number of Cores (--ncores):", self.expfact_ncores_spin)

        self.expfact_out_edit = QLineEdit("pfact_out")
        form_layout.addRow("Output File Prefix (--out):", self.expfact_out_edit)

        layout.addLayout(form_layout)
        run_button = QPushButton("Run ExPfact Fitting")
        run_button.setObjectName("run_button_expfact")
        run_button.clicked.connect(self._run_expfact_fitting)
        layout.addWidget(run_button)
        layout.addStretch()
        self.tab_widget.addTab(tab, "ExPfact Fitting")

    def _run_expfact_fitting(self):
        args = []
        if not self.expfact_dexp_edit.text() or \
           not self.expfact_ass_edit.text() or \
           not self.expfact_seq_edit.text():
            QMessageBox.warning(self, "Input Error", "Dexp, Assignments, and Sequence files are required.")
            return
        if not self.expfact_out_edit.text():
            QMessageBox.warning(self, "Input Error", "Output file prefix is required.")
            return

        args.extend(["--dexp", self.expfact_dexp_edit.text()])
        args.extend(["--ass", self.expfact_ass_edit.text()])
        args.extend(["--seq", self.expfact_seq_edit.text()])
        if self.expfact_pfact_cb.isChecked() and self.expfact_pfact_edit.text():
            args.extend(["--pfact", self.expfact_pfact_edit.text()])

        args.extend(["--temp", str(self.expfact_temp_spin.value())])
        args.extend(["--pH", str(self.expfact_ph_spin.value())])
        args.extend(["--harm", str(self.expfact_penalty_spin.value())])
        args.extend(["--tol", f"{self.expfact_tol_spin.value():.1e}"])
        if self.expfact_rand_cb.isChecked():
            args.extend(["--rand", str(self.expfact_rand_spin.value())])
        args.extend(["--rep", str(self.expfact_rep_spin.value())])
        args.extend(["--ncores", str(self.expfact_ncores_spin.value())])
        args.extend(["--out", self.expfact_out_edit.text()])
        
        self._run_script("exPfact.py", args)


    def _create_md2pfact_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        form_layout = QFormLayout()

        self.md2pfact_pdb_edit = QLineEdit()
        self.md2pfact_pdb_browse = QPushButton("Browse...")
        self.md2pfact_pdb_browse.clicked.connect(lambda: self._browse_file(self.md2pfact_pdb_edit, "PDB File", "PDB Files (*.pdb);;All Files (*.*)"))
        hb_pdb = QHBoxLayout()
        hb_pdb.addWidget(self.md2pfact_pdb_edit)
        hb_pdb.addWidget(self.md2pfact_pdb_browse)
        form_layout.addRow("PDB File (--pdb):", hb_pdb)

        self.md2pfact_dcd_edit = QLineEdit()
        self.md2pfact_dcd_browse = QPushButton("Browse...")
        self.md2pfact_dcd_browse.clicked.connect(lambda: self._browse_file(self.md2pfact_dcd_edit, "DCD File", "DCD Files (*.dcd);;All Files (*.*)"))
        hb_dcd = QHBoxLayout()
        hb_dcd.addWidget(self.md2pfact_dcd_edit)
        hb_dcd.addWidget(self.md2pfact_dcd_browse)
        form_layout.addRow("DCD File (--dcd):", hb_dcd)

        self.md2pfact_bc_spin = QDoubleSpinBox()
        self.md2pfact_bc_spin.setValue(0.35)
        form_layout.addRow("Beta_c (--bc):", self.md2pfact_bc_spin)

        self.md2pfact_bh_spin = QDoubleSpinBox()
        self.md2pfact_bh_spin.setValue(2.00)
        form_layout.addRow("Beta_h (--bh):", self.md2pfact_bh_spin)

        self.md2pfact_step_spin = QSpinBox()
        self.md2pfact_step_spin.setRange(1, 100000)
        self.md2pfact_step_spin.setValue(100)
        form_layout.addRow("Step (--step):", self.md2pfact_step_spin)

        self.md2pfact_out_edit = QLineEdit("md_pfact_results")
        form_layout.addRow("Output Prefix (--out):", self.md2pfact_out_edit)

        layout.addLayout(form_layout)
        run_button = QPushButton("Calculate Pfactors from MD")
        run_button.setObjectName("run_button_md2pfact")
        run_button.clicked.connect(self._run_md2pfact)
        layout.addWidget(run_button)

        self.md_plot_label = QLabel("Plot will appear here after successful run.")
        self.md_plot_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.md_plot_label.setMinimumSize(300,200)
        layout.addWidget(self.md_plot_label, 1)

        layout.addStretch(0)
        self.tab_widget.addTab(tab, "MD to Pfact")

    def _run_md2pfact(self):
        args = []
        if not self.md2pfact_pdb_edit.text() or not self.md2pfact_dcd_edit.text():
            QMessageBox.warning(self, "Input Error", "PDB and DCD files are required.")
            return
        if not self.md2pfact_out_edit.text():
            QMessageBox.warning(self, "Input Error", "Output Prefix is required.")
            return
        
        args.extend(["--pdb", self.md2pfact_pdb_edit.text()])
        args.extend(["--dcd", self.md2pfact_dcd_edit.text()])
        args.extend(["--out", self.md2pfact_out_edit.text()])
        args.extend(["--bc", str(self.md2pfact_bc_spin.value())])
        args.extend(["--bh", str(self.md2pfact_bh_spin.value())])
        args.extend(["--step", str(self.md2pfact_step_spin.value())])

        self.md_plot_label.setText("Plot will appear here after successful run.")
        self.md_plot_label.setPixmap(QPixmap())

        self._run_script("MD2Pfact.py", args)

    def _display_md_plot(self, image_path):
        pixmap = QPixmap(image_path)
        if pixmap.isNull():
            self.md_plot_label.setText(f"Failed to load plot: {image_path}")
            self._log_message(f"Failed to load plot: {image_path}", is_error=True)
        else:
            scaled_pixmap = pixmap.scaled(self.md_plot_label.size(),
                                          Qt.AspectRatioMode.KeepAspectRatio,
                                          Qt.TransformationMode.SmoothTransformation)
            self.md_plot_label.setPixmap(scaled_pixmap)
            self._log_message(f"Displayed plot: {image_path}")


    def _create_pfact2dpred_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        form_layout = QFormLayout()

        self.pfact2dpred_pfact_edit = QLineEdit()
        self.pfact2dpred_pfact_browse = QPushButton("Browse...")
        self.pfact2dpred_pfact_browse.clicked.connect(lambda: self._browse_file(self.pfact2dpred_pfact_edit, "Pfact File", "Pfact Files (*.pfact *.txt);;All Files (*.*)"))
        hb_pfact = QHBoxLayout()
        hb_pfact.addWidget(self.pfact2dpred_pfact_edit)
        hb_pfact.addWidget(self.pfact2dpred_pfact_browse)
        form_layout.addRow("Pfact File (--pfact):", hb_pfact)

        self.pfact2dpred_ass_edit = QLineEdit()
        self.pfact2dpred_ass_browse = QPushButton("Browse...")
        self.pfact2dpred_ass_browse.clicked.connect(lambda: self._browse_file(self.pfact2dpred_ass_edit, "Assignments File", "List Files (*.ass *.list *.txt);;All Files (*.*)"))
        hb_ass = QHBoxLayout()
        hb_ass.addWidget(self.pfact2dpred_ass_edit)
        hb_ass.addWidget(self.pfact2dpred_ass_browse)
        form_layout.addRow("Assignments File (--ass):", hb_ass)

        self.pfact2dpred_times_edit = QLineEdit()
        self.pfact2dpred_times_browse = QPushButton("Browse...")
        self.pfact2dpred_times_browse.clicked.connect(lambda: self._browse_file(self.pfact2dpred_times_edit, "Times File", "Times Files (*.times *.txt);;All Files (*.*)"))
        hb_times = QHBoxLayout()
        hb_times.addWidget(self.pfact2dpred_times_edit)
        hb_times.addWidget(self.pfact2dpred_times_browse)
        form_layout.addRow("Times File (--times):", hb_times)

        self.pfact2dpred_seq_edit = QLineEdit()
        self.pfact2dpred_seq_browse = QPushButton("Browse...")
        self.pfact2dpred_seq_browse.clicked.connect(lambda: self._browse_file(self.pfact2dpred_seq_edit, "Sequence File", "Sequence Files (*.seq *.fasta *.txt);;All Files (*.*)"))
        hb_seq = QHBoxLayout()
        hb_seq.addWidget(self.pfact2dpred_seq_edit)
        hb_seq.addWidget(self.pfact2dpred_seq_browse)
        form_layout.addRow("Sequence File (--seq):", hb_seq)

        self.pfact2dpred_temp_spin = QDoubleSpinBox()
        self.pfact2dpred_temp_spin.setRange(0, 500); self.pfact2dpred_temp_spin.setValue(300)
        form_layout.addRow("Temperature (--temp, K):", self.pfact2dpred_temp_spin)

        self.pfact2dpred_ph_spin = QDoubleSpinBox()
        self.pfact2dpred_ph_spin.setRange(0, 14); self.pfact2dpred_ph_spin.setValue(7.0)
        form_layout.addRow("pH (--pH):", self.pfact2dpred_ph_spin)

        self.pfact2dpred_nrep_spin = QSpinBox()
        self.pfact2dpred_nrep_spin.setRange(1, 100); self.pfact2dpred_nrep_spin.setValue(1)
        form_layout.addRow("Num Replicates (--nrep):", self.pfact2dpred_nrep_spin)

        self.pfact2dpred_eps_spin = QDoubleSpinBox()
        self.pfact2dpred_eps_spin.setRange(0, 1); self.pfact2dpred_eps_spin.setDecimals(4)
        self.pfact2dpred_eps_spin.setValue(0.0)
        form_layout.addRow("Epsilon (--eps):", self.pfact2dpred_eps_spin)

        self.pfact2dpred_out_edit = QLineEdit("dpred_out")
        form_layout.addRow("Output File Prefix (--out):", self.pfact2dpred_out_edit)

        layout.addLayout(form_layout)
        run_button = QPushButton("Predict D-uptake")
        run_button.setObjectName("run_button_pfact2dpred")
        run_button.clicked.connect(self._run_pfact2dpred)
        layout.addWidget(run_button)
        layout.addStretch()
        self.tab_widget.addTab(tab, "Pfact to Dpred")

    def _run_pfact2dpred(self):
        args = []
        if not all([self.pfact2dpred_pfact_edit.text(), self.pfact2dpred_ass_edit.text(),
                    self.pfact2dpred_times_edit.text(), self.pfact2dpred_seq_edit.text()]):
            QMessageBox.warning(self, "Input Error", "Pfact, Assignments, Times, and Sequence files are required.")
            return
        if not self.pfact2dpred_out_edit.text():
            QMessageBox.warning(self, "Input Error", "Output file prefix is required."); return

        args.extend(["--pfact", self.pfact2dpred_pfact_edit.text()])
        args.extend(["--ass", self.pfact2dpred_ass_edit.text()])
        args.extend(["--times", self.pfact2dpred_times_edit.text()])
        args.extend(["--seq", self.pfact2dpred_seq_edit.text()])
        args.extend(["--out", self.pfact2dpred_out_edit.text()])
        args.extend(["--temp", str(self.pfact2dpred_temp_spin.value())])
        args.extend(["--pH", str(self.pfact2dpred_ph_spin.value())])
        args.extend(["--nrep", str(self.pfact2dpred_nrep_spin.value())])
        args.extend(["--eps", str(self.pfact2dpred_eps_spin.value())])
        
        self._run_script("pfact2dpred.py", args)


    def _create_isotopes_tab(self):
        tab = QWidget()
        main_isotopes_layout = QVBoxLayout(tab)

        hisotope_group = QWidget()
        hisotope_layout = QVBoxLayout(hisotope_group)
        hisotope_layout.addWidget(QLabel("<b>Fully Protonated Envelope (Hisotope.py)</b>"))
        hisotope_form = QFormLayout()

        self.hisotope_seq_edit = QLineEdit("SAMPLE")
        hisotope_form.addRow("Sequence (--seq):", self.hisotope_seq_edit)

        self.hisotope_z_spin = QSpinBox()
        self.hisotope_z_spin.setRange(1, 10); self.hisotope_z_spin.setValue(3)
        hisotope_form.addRow("Charge State (--z):", self.hisotope_z_spin)

        hisotope_layout.addLayout(hisotope_form)
        run_hisotope_button = QPushButton("Calculate Fully Protonated Envelope")
        run_hisotope_button.setObjectName("run_button_hisotope")
        run_hisotope_button.clicked.connect(self._run_hisotope)
        hisotope_layout.addWidget(run_hisotope_button)
        main_isotopes_layout.addWidget(hisotope_group)

        isenv_group = QWidget()
        isenv_layout = QVBoxLayout(isenv_group)
        isenv_layout.addWidget(QLabel("<b>Time-Dependent Isotopic Envelope (isotopic_envelope.py)</b>"))
        isenv_form = QFormLayout()

        self.isenv_mode_combo = QComboBox()
        self.isenv_mode_combo.addItems(["Predict (p)", "Compare (c)"])
        isenv_form.addRow("Mode (--mode):", self.isenv_mode_combo)

        self.isenv_ass_edit = QLineEdit()
        self.isenv_ass_browse = QPushButton("Browse...")
        self.isenv_ass_browse.clicked.connect(lambda: self._browse_file(self.isenv_ass_edit, "Assignments File", "List Files (*.ass *.list *.txt);;All Files (*.*)"))
        hb_ass = QHBoxLayout()
        hb_ass.addWidget(self.isenv_ass_edit)
        hb_ass.addWidget(self.isenv_ass_browse)
        isenv_form.addRow("Assignments File (--ass):", hb_ass)

        self.isenv_seq_edit = QLineEdit()
        self.isenv_seq_browse = QPushButton("Browse...")
        self.isenv_seq_browse.clicked.connect(lambda: self._browse_file(self.isenv_seq_edit, "Sequence File", "Sequence Files (*.seq *.fasta *.txt);;All Files (*.*)"))
        hb_seq = QHBoxLayout()
        hb_seq.addWidget(self.isenv_seq_edit)
        hb_seq.addWidget(self.isenv_seq_browse)
        isenv_form.addRow("Sequence File (--seq):", hb_seq)

        self.isenv_pfact_edit = QLineEdit()
        self.isenv_pfact_browse = QPushButton("Browse...")
        self.isenv_pfact_browse.clicked.connect(lambda: self._browse_file(self.isenv_pfact_edit, "Pfact File", "Pfact Files (*.pfact *.txt);;All Files (*.*)"))
        hb_pfact = QHBoxLayout()
        hb_pfact.addWidget(self.isenv_pfact_edit)
        hb_pfact.addWidget(self.isenv_pfact_browse)
        isenv_form.addRow("Pfact File (--pfact):", hb_pfact)

        self.isenv_times_edit = QLineEdit()
        self.isenv_times_browse = QPushButton("Browse...")
        self.isenv_times_browse.clicked.connect(lambda: self._browse_file(self.isenv_times_edit, "Times File", "Times Files (*.times *.txt);;All Files (*.*)"))
        hb_times = QHBoxLayout()
        hb_times.addWidget(self.isenv_times_edit)
        hb_times.addWidget(self.isenv_times_browse)
        isenv_form.addRow("Times File (--times):", hb_times)

        self.isenv_T_label_spin = QDoubleSpinBox(); self.isenv_T_label_spin.setValue(300)
        isenv_form.addRow("Labeling Temp (--T_label, K):", self.isenv_T_label_spin)
        self.isenv_pH_label_spin = QDoubleSpinBox(); self.isenv_pH_label_spin.setValue(7.0)
        isenv_form.addRow("Labeling pH (--pH_label):", self.isenv_pH_label_spin)

        self.isenv_pep_spin = QSpinBox(); self.isenv_pep_spin.setRange(1,10000); self.isenv_pep_spin.setValue(1)
        isenv_form.addRow("Peptide Index (--pep):", self.isenv_pep_spin)
        self.isenv_z_spin = QSpinBox(); self.isenv_z_spin.setRange(1,10); self.isenv_z_spin.setValue(3)
        isenv_form.addRow("Charge State (--z):", self.isenv_z_spin)

        self.isenv_prefix_edit = QLineEdit("isotope_env_results")
        isenv_form.addRow("Output Prefix (--prefix):", self.isenv_prefix_edit)

        self.isenv_T_quench_label = QLabel("Quench Temp (--T_quench, K):")
        self.isenv_T_quench_spin = QDoubleSpinBox(); self.isenv_T_quench_spin.setValue(277)
        self.isenv_pH_quench_label = QLabel("Quench pH (--pH_quench):")
        self.isenv_pH_quench_spin = QDoubleSpinBox(); self.isenv_pH_quench_spin.setValue(2.5)

        isenv_form.addRow(self.isenv_T_quench_label, self.isenv_T_quench_spin)
        isenv_form.addRow(self.isenv_pH_quench_label, self.isenv_pH_quench_spin)

        def _toggle_compare_fields(index):
            is_compare = self.isenv_mode_combo.currentText().startswith("Compare")
            self.isenv_T_quench_label.setVisible(is_compare)
            self.isenv_T_quench_spin.setVisible(is_compare)
            self.isenv_pH_quench_label.setVisible(is_compare)
            self.isenv_pH_quench_spin.setVisible(is_compare)
        self.isenv_mode_combo.currentIndexChanged.connect(_toggle_compare_fields)
        _toggle_compare_fields(0)

        isenv_layout.addLayout(isenv_form)
        run_isenv_button = QPushButton("Calculate/Compare Envelopes")
        run_isenv_button.setObjectName("run_button_isenv")
        run_isenv_button.clicked.connect(self._run_isotopic_envelope)
        isenv_layout.addWidget(run_isenv_button)

        main_isotopes_layout.addWidget(isenv_group)
        main_isotopes_layout.addStretch()
        self.tab_widget.addTab(tab, "Isotopic Envelopes")

    def _run_hisotope(self):
        args = []
        if not self.hisotope_seq_edit.text():
            QMessageBox.warning(self, "Input Error", "Sequence is required for Hisotope."); return
        args.extend(["--seq", self.hisotope_seq_edit.text()])
        args.extend(["--z", str(self.hisotope_z_spin.value())])
        self._run_script("Hisotope.py", args)

    def _run_isotopic_envelope(self):
        args = []
        if not all([self.isenv_ass_edit.text(), self.isenv_seq_edit.text(),
                    self.isenv_pfact_edit.text(), self.isenv_times_edit.text()]):
            QMessageBox.warning(self, "Input Error", "Ass, Seq, Pfact, and Times files are required.")
            return
        if not self.isenv_prefix_edit.text():
            QMessageBox.warning(self, "Input Error", "Output Prefix is required."); return

        mode = "p" if self.isenv_mode_combo.currentText().startswith("Predict") else "c"
        args.extend(["--mode", mode])
        args.extend(["--ass", self.isenv_ass_edit.text()])
        args.extend(["--seq", self.isenv_seq_edit.text()])
        args.extend(["--pfact", self.isenv_pfact_edit.text()])
        args.extend(["--times", self.isenv_times_edit.text()])
        args.extend(["--prefix", self.isenv_prefix_edit.text()])
        args.extend(["--T_label", str(self.isenv_T_label_spin.value())])
        args.extend(["--pH_label", str(self.isenv_pH_label_spin.value())])
        args.extend(["--pep", str(self.isenv_pep_spin.value())])
        args.extend(["--z", str(self.isenv_z_spin.value())])

        if mode == 'c':
            args.extend(["--T_quench", str(self.isenv_T_quench_spin.value())])
            args.extend(["--pH_quench", str(self.isenv_pH_quench_spin.value())])
        
        self._run_script("isotopic_envelope.py", args)


    def _create_utilities_tab(self):
        utilities_main_tab = QWidget()
        util_tab_widget = QTabWidget(utilities_main_tab)
        main_layout = QVBoxLayout(utilities_main_tab)
        main_layout.addWidget(util_tab_widget)

        # --- Sub-Tab 1: Process DynamX Cluster ---
        dnyx_tab = QWidget()
        dnyx_layout = QVBoxLayout(dnyx_tab)
        dnyx_form = QFormLayout()

        self.dnyx_csv_edit = QLineEdit()
        self.dnyx_csv_browse = QPushButton("Browse...")
        self.dnyx_csv_browse.clicked.connect(lambda: self._browse_file(self.dnyx_csv_edit, "DynamX CSV File", "CSV Files (*.csv);;All Files (*.*)"))
        hb_csv = QHBoxLayout()
        hb_csv.addWidget(self.dnyx_csv_edit)
        hb_csv.addWidget(self.dnyx_csv_browse)
        dnyx_form.addRow("DynamX .csv File:", hb_csv)

        self.dnyx_norm_method_combo = QComboBox()
        self.dnyx_norm_method_combo.addItems(["1: Theoretical MaxUptake", "2: Fully Deuterated Sample"])
        dnyx_form.addRow("Normalization Method:", self.dnyx_norm_method_combo)

        self.dnyx_d2o_perc_label = QLabel("D2O Percentage (0-100):")
        self.dnyx_d2o_perc_spin = QSpinBox()
        self.dnyx_d2o_perc_spin.setRange(0, 100); self.dnyx_d2o_perc_spin.setValue(90)
        dnyx_form.addRow(self.dnyx_d2o_perc_label, self.dnyx_d2o_perc_spin)

        def _toggle_d2o_perc(index):
            self.dnyx_d2o_perc_label.setVisible(index == 0)
            self.dnyx_d2o_perc_spin.setVisible(index == 0)
        self.dnyx_norm_method_combo.currentIndexChanged.connect(_toggle_d2o_perc)
        _toggle_d2o_perc(0)

        dnyx_layout.addLayout(dnyx_form)
        run_dnyx_button = QPushButton("Process DynamX File")
        run_dnyx_button.setObjectName("run_button_dnyx")
        run_dnyx_button.clicked.connect(self._run_process_dnyx)
        dnyx_layout.addWidget(run_dnyx_button)
        dnyx_layout.addStretch()
        util_tab_widget.addTab(dnyx_tab, "Process DynamX")

        # --- Sub-Tab 2: Cross Validation ---
        cv_tab = QWidget()
        cv_layout = QVBoxLayout(cv_tab)
        cv_form = QFormLayout()

        self.cv_dexp_edit = QLineEdit()
        self.cv_dexp_browse = QPushButton("Browse...")
        self.cv_dexp_browse.clicked.connect(lambda: self._browse_file(self.cv_dexp_edit, "Dexp File", "Data Files (*.dexp *.txt *.dat);;All Files (*.*)"))
        hb_cv_dexp = QHBoxLayout()
        hb_cv_dexp.addWidget(self.cv_dexp_edit)
        hb_cv_dexp.addWidget(self.cv_dexp_browse)
        cv_form.addRow("Dexp File (--dexp):", hb_cv_dexp)

        self.cv_ass_edit = QLineEdit()
        self.cv_ass_browse = QPushButton("Browse...")
        self.cv_ass_browse.clicked.connect(lambda: self._browse_file(self.cv_ass_edit, "Assignments File", "List Files (*.ass *.list *.txt);;All Files (*.*)"))
        hb_cv_ass = QHBoxLayout()
        hb_cv_ass.addWidget(self.cv_ass_edit)
        hb_cv_ass.addWidget(self.cv_ass_browse)
        cv_form.addRow("Assignments File (--ass):", hb_cv_ass)

        self.cv_seq_edit = QLineEdit()
        self.cv_seq_browse = QPushButton("Browse...")
        self.cv_seq_browse.clicked.connect(lambda: self._browse_file(self.cv_seq_edit, "Sequence File", "Sequence Files (*.seq *.fasta *.txt);;All Files (*.*)"))
        hb_cv_seq = QHBoxLayout()
        hb_cv_seq.addWidget(self.cv_seq_edit)
        hb_cv_seq.addWidget(self.cv_seq_browse)
        cv_form.addRow("Sequence File (--seq):", hb_cv_seq)

        self.cv_temp_spin = QDoubleSpinBox(); self.cv_temp_spin.setRange(0,500); self.cv_temp_spin.setValue(300)
        cv_form.addRow("Temperature (--temp, K):", self.cv_temp_spin)
        self.cv_ph_spin = QDoubleSpinBox(); self.cv_ph_spin.setRange(0,14); self.cv_ph_spin.setValue(7.0)
        cv_form.addRow("pH (--pH):", self.cv_ph_spin)

        cv_layout.addLayout(cv_form)
        run_cv_button = QPushButton("Run Cross Validation")
        run_cv_button.setObjectName("run_button_cv")
        run_cv_button.clicked.connect(self._run_cross_validation)
        cv_layout.addWidget(run_cv_button)
        cv_layout.addStretch()
        util_tab_widget.addTab(cv_tab, "Cross Validation")

        # --- Sub-Tab 3: Descriptive Statistics ---
        desc_stats_tab = QWidget()
        desc_stats_layout = QVBoxLayout(desc_stats_tab)
        desc_stats_form = QFormLayout()

        self.desc_stats_res_edit = QLineEdit("pfact_out")
        desc_stats_form.addRow("Results Prefix (--res):", self.desc_stats_res_edit)

        self.desc_stats_top_spin = QSpinBox()
        self.desc_stats_top_spin.setRange(1,100); self.desc_stats_top_spin.setValue(50)
        desc_stats_form.addRow("Top % (--top):", self.desc_stats_top_spin)

        desc_stats_layout.addLayout(desc_stats_form)
        run_desc_stats_button = QPushButton("Calculate Descriptive Statistics")
        run_desc_stats_button.setObjectName("run_button_desc_stats")
        run_desc_stats_button.clicked.connect(self._run_descriptive_stats)
        desc_stats_layout.addWidget(run_desc_stats_button)
        desc_stats_layout.addStretch()
        util_tab_widget.addTab(desc_stats_tab, "Descriptive Statistics")

        # --- Sub-Tab 4: Clustering  ---
        cluster_tab = QWidget()
        cluster_layout = QVBoxLayout(cluster_tab)
        cluster_form = QFormLayout()
        
        self.cluster_ass_edit = QLineEdit()
        self.cluster_ass_browse = QPushButton("Browse...")
        self.cluster_ass_browse.clicked.connect(lambda: self._browse_file(
            self.cluster_ass_edit, "Assignments File", "Assignment Files (*.ass *.list);;All Files (*.*)"
        ))
        hb_cluster_ass = QHBoxLayout()
        hb_cluster_ass.addWidget(self.cluster_ass_edit)
        hb_cluster_ass.addWidget(self.cluster_ass_browse)
        cluster_form.addRow("Assignments File (--ass):", hb_cluster_ass)

        # Using self.cluster_all_sp_edit to match clustering.py --all_sp argument
        self.cluster_all_sp_edit = QLineEdit() 
        self.cluster_all_sp_browse = QPushButton("Browse...")
        self.cluster_all_sp_browse.clicked.connect(lambda: self._browse_file(
            self.cluster_all_sp_edit, "all.sp Data File (Optional)", "SP Files (*.sp);;Data Files (*.dat *.txt);;All Files (*.*)"
        ))
        hb_cluster_all_sp = QHBoxLayout()
        hb_cluster_all_sp.addWidget(self.cluster_all_sp_edit)
        hb_cluster_all_sp.addWidget(self.cluster_all_sp_browse)
        # Label reflects the argument name used in clustering.py
        cluster_form.addRow("all.sp Data File (--all_sp):", hb_cluster_all_sp) 
        
        cluster_layout.addLayout(cluster_form)
        run_cluster_button = QPushButton("Run Clustering (requires R & mclust)")
        run_cluster_button.setObjectName("run_button_clustering")
        run_cluster_button.clicked.connect(self._run_clustering)
        cluster_layout.addWidget(run_cluster_button)
        cluster_layout.addStretch()
        util_tab_widget.addTab(cluster_tab, "Clustering")

        self.tab_widget.addTab(utilities_main_tab, "Utilities")


    def _run_process_dnyx(self):
        args = []
        if not self.dnyx_csv_edit.text():
            QMessageBox.warning(self, "Input Error", "DynamX CSV file is required.")
            return
        
        args.append(self.dnyx_csv_edit.text())
        args.extend(["--norm_method", self.dnyx_norm_method_combo.currentText().split(':')[0]])
        if self.dnyx_norm_method_combo.currentIndex() == 0: 
             args.extend(["--d2o_perc", str(self.dnyx_d2o_perc_spin.value())])
        
        self._run_script("process_DnXcluster.py", args)

    def _run_cross_validation(self):
        args = []
        if not all([self.cv_dexp_edit.text(), self.cv_ass_edit.text(), self.cv_seq_edit.text()]):
             QMessageBox.warning(self, "Input Error", "Dexp, Ass, and Seq files are required for CV.")
             return
        
        args.extend(["--dexp", self.cv_dexp_edit.text()])
        args.extend(["--ass", self.cv_ass_edit.text()])
        args.extend(["--seq", self.cv_seq_edit.text()])
        args.extend(["--temp", str(self.cv_temp_spin.value())])
        args.extend(["--pH", str(self.cv_ph_spin.value())])
        
        self._run_script("cross_validation.py", args)

    def _run_descriptive_stats(self):
        args = []
        if not self.desc_stats_res_edit.text():
            QMessageBox.warning(self, "Input Error", "Results Prefix is required.")
            return
        
        args.extend(["--res", self.desc_stats_res_edit.text()])
        args.extend(["--top", str(self.desc_stats_top_spin.value())])
        
        self._run_script("descriptive_statistics.py", args)

    def _run_clustering(self):
        args = []
        
        if not self.cluster_ass_edit.text():
            QMessageBox.warning(self, "Input Error", "Assignments file (--ass) is required for clustering.")
            return
        args.extend(["--ass", self.cluster_ass_edit.text()])

        if hasattr(self, 'cluster_all_sp_edit') and self.cluster_all_sp_edit.text():
            all_sp_path = self.cluster_all_sp_edit.text()
            if not os.path.isfile(all_sp_path):
                QMessageBox.warning(self, "Input Error",
                                    f"The specified 'all.sp' file was not found: {all_sp_path}\n"
                                    "Clustering will proceed assuming 'all.sp' is in the 'results' directory or not needed if path is handled differently by the R script.")
            args.extend(["--all_sp", all_sp_path]) # Use --all_sp
        
        log_msg = "Ensure R and mclust are correctly set up. "
        if not (hasattr(self, 'cluster_all_sp_edit') and self.cluster_all_sp_edit.text()):
            log_msg += "Since 'all.sp' path is not specified, 'clustering.py' will expect a file named 'all.sp' in the 'results' directory if needed by its R script logic."
        
        self._log_message(log_msg)
        self._run_script("clustering.py", args)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    
    if not os.path.exists(SCRIPT_DIR):
        print(f"Warning: Python script directory '{SCRIPT_DIR}' not found. Creating it.")
        try:
            os.makedirs(SCRIPT_DIR)
            print(f"Created directory '{SCRIPT_DIR}'. Place your ExPfact Python scripts there.")
        except OSError as e:
            print(f"Critical Error: Could not create script directory {SCRIPT_DIR}: {e}. Exiting.")
            sys.exit(1)
    else:
        print(f"Using Python scripts from '{SCRIPT_DIR}'.")

    expected_r_dir_path = os.path.abspath(os.path.join(SCRIPT_DIR, R_SCRIPT_DIR_RELATIVE_TO_PY_SCRIPT))
    if not os.path.exists(expected_r_dir_path):
         print(f"Warning: R script directory for clustering not found at expected location: {expected_r_dir_path}")
         print(f"Clustering functionality might fail if 'multi.r' is not in '{expected_r_dir_path}'.")
    else:
        print(f"R scripts for clustering expected in '{expected_r_dir_path}'.")

    main_win = ExPfactGUI()
    main_win.show()
    sys.exit(app.exec())
