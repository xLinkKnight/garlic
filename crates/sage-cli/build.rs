//! Embeds the XNova icon in `sage.exe` (the same `assets/xnova.ico` as
//! `sage-gui.exe`), so Explorer and shortcuts show it. A missing resource
//! compiler only costs the icon.

const ICON: &str = "../sage-gui/assets/xnova.ico";

fn main() {
    println!("cargo:rerun-if-changed={ICON}");
    #[cfg(windows)]
    {
        let mut res = winresource::WindowsResource::new();
        res.set_icon(ICON);
        if let Err(e) = res.compile() {
            println!("cargo:warning=sage.exe built without its icon: {e}");
        }
    }
}
