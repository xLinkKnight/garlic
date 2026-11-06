//! Embeds the application icon (`assets/xnova.ico`) in `sage-gui.exe`, so
//! Explorer, shortcuts and the taskbar show it. A missing resource compiler
//! only costs the file icon: the window icon is set at run time anyway.

fn main() {
    println!("cargo:rerun-if-changed=assets/xnova.ico");
    #[cfg(windows)]
    {
        let mut res = winresource::WindowsResource::new();
        res.set_icon("assets/xnova.ico");
        if let Err(e) = res.compile() {
            println!("cargo:warning=sage-gui.exe built without its icon: {e}");
        }
    }
}
