/*
 * Privacy-conscious analytics for scoring.hellstrom.pw.
 *
 * Setup:
 *   1. Create a GoatCounter site at https://www.goatcounter.com/
 *   2. Replace the empty GOATCOUNTER_CODE below with your site code.
 *      Example: const GOATCOUNTER_CODE = "scoring";
 *
 * Until a code is configured this file exits without making any network
 * requests, so it is safe to deploy before the analytics account exists.
 */
(() => {
  "use strict";

  const GOATCOUNTER_CODE = "";
  const SHOW_VISITOR_COUNT = true;

  if (!GOATCOUNTER_CODE) return;

  const endpoint = `https://${GOATCOUNTER_CODE}.goatcounter.com/count`;

  // Keep campaign parameters available to GoatCounter while avoiding arbitrary
  // query strings in page names. The canonical URL already identifies the page.
  window.goatcounter = {
    path: () => location.pathname || "/",
    referrer: () => {
      const params = new URLSearchParams(location.search);
      return (
        params.get("ref") ||
        params.get("utm_campaign") ||
        params.get("utm_source") ||
        document.referrer ||
        undefined
      );
    },
  };

  const tracker = document.createElement("script");
  tracker.async = true;
  tracker.src = "https://gc.zgo.at/count.js";
  tracker.dataset.goatcounter = endpoint;

  tracker.addEventListener("load", () => {
    // Record file downloads and outbound clicks as events in addition to the
    // normal page-view/referrer data.
    document.querySelectorAll("a[href]").forEach((link) => {
      const href = link.getAttribute("href");
      if (!href || href.startsWith("#") || href.startsWith("mailto:") || href.startsWith("tel:")) return;

      let url;
      try {
        url = new URL(href, location.href);
      } catch (_) {
        return;
      }

      const extensionMatch = url.pathname.match(/\.([a-z0-9]{1,8})$/i);
      const extension = extensionMatch ? extensionMatch[1].toLowerCase() : "";
      const downloadTypes = new Set(["pdf", "csv", "r", "odt", "docx", "zip", "png", "jpg", "jpeg", "txt"]);
      const isDownload = url.origin === location.origin && downloadTypes.has(extension);
      const isOutbound = url.origin !== location.origin;

      if (!isDownload && !isOutbound) return;

      const eventName = isDownload
        ? `download-${url.pathname.replace(/^\/+/, "")}`
        : `outbound-${url.hostname}`;

      link.addEventListener("click", () => {
        if (!window.goatcounter || typeof window.goatcounter.count !== "function") return;
        window.goatcounter.count({
          path: eventName,
          title: (link.textContent || link.getAttribute("aria-label") || href).trim().slice(0, 200),
          event: true,
          no_session: true,
        });
      });
    });

    if (!SHOW_VISITOR_COUNT || !window.goatcounter || typeof window.goatcounter.visit_count !== "function") return;

    const footer = document.querySelector("footer");
    if (!footer) return;

    const counter = document.createElement("span");
    counter.id = "visitor-count";
    counter.setAttribute("aria-label", "Page views");
    counter.style.marginLeft = ".75rem";
    footer.appendChild(counter);

    // GoatCounter requires the site's "Allow adding visitor counts" setting
    // to be enabled. If it is disabled the normal analytics still works.
    try {
      window.goatcounter.visit_count({ append: "#visitor-count", type: "text" });
    } catch (_) {
      counter.remove();
    }
  });

  document.head.appendChild(tracker);
})();
